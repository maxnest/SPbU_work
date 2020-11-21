try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="File with phylostratigraphic levels described")
parser.add_argument('--go_univ_merged', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of annotation with eggNOG (GeneOntologyUniverse for topGO)."
                         "The last column should contain information about species")
parser.add_argument('--phylostratr_merged', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of phylostratr analysis."
                         "The last column should contain information about species")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def levels_parsing(levels, levels_dict):
    for line in levels:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def go_univ_parsing(go_univ_merged, go_univ_dict):
    header = go_univ_merged.readline()
    for line in go_univ_merged:
        description = line.strip().split("\t")
        protein_ID, terms, species = description[0], description[1], description[2]
        if species not in go_univ_dict.keys():
            go_univ_dict[species] = []
        go_univ_dict[species].append(protein_ID)


def phylostratr_parsing(phylostratr_merged, phylostratr_dict):
    header = phylostratr_merged.readline()
    for line in phylostratr_merged:
        description = line.strip().split("\t")
        protein_ID, phylostrata, species = description[0], description[-2][1:-1], description[-1]
        if species not in phylostratr_dict.keys():
            phylostratr_dict[species] = {}
        phylostratr_dict[species][protein_ID] = phylostrata


def generalization(levels_dict, go_univ_dict, phylostratr_dict, summary_dict):
    species_list = [key for key in go_univ_dict.keys()]

    for phylostrata, level in levels_dict.items():
        summary_dict[level] = {species: [] for species in species_list}

    for species, protein_list in go_univ_dict.items():
        for protein in protein_list:
            summary_dict[levels_dict[phylostratr_dict[species][protein]]][species].append(protein)


def output_writing(output, summary_dict, go_univ_dict):
    species_list = [key for key in go_univ_dict.keys()]

    with open("{output}_proteins_with_GO-terms_in_phylostrates.tsv".format(output=output), 'a') as output_file:
        output_file.write("Levels\t{species}\n".format(species="\t".join(species_list)))
        for level, species_values in summary_dict.items():
            values = ["{len} ({percent}%)".format(len=len(species_values[species]),
                      percent=round((len(species_values[species])/len(go_univ_dict[species]))*100, 2))
                      for species in species_list]
            output_file.write("{level}\t{values}\n".format(level=level, values="\t".join(values)))


if __name__ == "__main__":
    levels_dict, go_univ_dict, phylostratr_dict, summary_dict = {}, {}, {}, {}
    levels_parsing(args.levels, levels_dict)
    go_univ_parsing(args.go_univ_merged, go_univ_dict)
    phylostratr_parsing(args.phylostratr_merged, phylostratr_dict)
    generalization(levels_dict, go_univ_dict, phylostratr_dict, summary_dict)
    output_writing(args.output, summary_dict, go_univ_dict)





