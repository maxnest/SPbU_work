try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--matrix', type=argparse.FileType('r'), required=True,
                    help="Original matrix with results of primary metabolism summary")
parser.add_argument('--apro_eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with annotation results created by eggNOG-mapper for Amoeboaphelidium protococcarum "
                         "proteins")
parser.add_argument('--ains_eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with annotation results created by eggNOG-mapper for Aphelidium insulamus "
                         "proteins")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def original_table_parsing(matrix, matrix_dict, species_list):
    species_list.extend(matrix.readline().strip().split("\t")[8:])
    print("Species TAGs in the input matrix: {species}".format(species=" ".join(species_list)))

    for line in matrix:
        description = line.strip().split("\t")
        ko, values = description[4], description[8:]
        matrix_dict[ko] = {species: 0 for species in species_list}
        for species in species_list:
            matrix_dict[ko][species] += int(values[species_list.index(species)])


def annotation_parsing(annot, annot_dict):
    for line in annot:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            pep_ID, values = description[0], description[1:]
            if len(values[7]) != 0: # KEGG_ko
                for ko in values[7].split(","):
                    ko_ID = ko.split(":")[1]
                    if ko_ID not in annot_dict.keys():
                        annot_dict[ko_ID] = []
                    annot_dict[ko_ID].append(pep_ID)


def append_counts_in_matrix(matrix_dict, annot_dict, species_tag):
    for ko, values in matrix_dict.items():
        if ko in annot_dict.keys():
            values[species_tag] = len(set(annot_dict[ko]))
        else:
            values[species_tag] = 0


def output_matrix_creating(species_list, matrix_dict, apro_annot_dict, ains_annot_dict, output):
    append_counts_in_matrix(matrix_dict, apro_annot_dict, "Amo_pr")
    append_counts_in_matrix(matrix_dict, ains_annot_dict, "Aph_in")
    species_list.extend(["Amo_pr", "Aph_in"])

    with open("{output}.tsv".format(output=output), 'a') as output_matrix:
        output_matrix.write("KOs\t{species}\n".format(species="\t".join(species_list)))
        for ko, values in matrix_dict.items():
            species_values = [str(values[species]) for species in species_list]
            output_matrix.write("{ko}\t{species_values}\n".format(ko=ko,
                                                                  species_values="\t".join(species_values)))


if __name__ == "__main__":
    species_list, matrix_dict, apro_annot_dict, ains_annot_dict = [], {}, {}, {}
    print("***** Input files parsing *****")
    original_table_parsing(args.matrix, matrix_dict, species_list)
    annotation_parsing(args.apro_eggnog, apro_annot_dict)
    annotation_parsing(args.ains_eggnog, ains_annot_dict)
    print("***** Output file creating *****")
    output_matrix_creating(species_list, matrix_dict, apro_annot_dict, ains_annot_dict, args.output)
