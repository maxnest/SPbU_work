try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="File with phylostratigraphic levels described")
parser.add_argument('--exsec_merged', type=argparse.FileType('r'), required=True,
                    help="Merged table with potential excretory/secretory sequences IDs."
                         "The last column should contain information about species")
parser.add_argument('--phylostratr_merged', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of phylostratr analysis."
                         "The last column should contain information about species")
parser.add_argument('--exsec_tag', type=str, required=True,
                    help="Tag for set of excretory/secretory sequences")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def levels_parsing(levels, levels_dict):
    for line in levels:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def phylostratr_parsing(phylostratr_merged, phylostratr_dict):
    header = phylostratr_merged.readline()
    for line in phylostratr_merged:
        description = line.strip().split("\t")
        protein_ID, phylostrata, species = description[0], description[-2][1:-1], description[-1]
        if species not in phylostratr_dict.keys():
            phylostratr_dict[species] = {}
        phylostratr_dict[species][protein_ID] = phylostrata


def exsec_parsing(exsec_merged, exsec_dict):
    header = exsec_merged.readline()
    for line in exsec_merged:
        description = line.strip().split("\t")
        protein_ID, species = description[0], description[-1]
        if species not in exsec_dict.keys():
            exsec_dict[species] = []
        exsec_dict[species].append(protein_ID)


def generalization(levels_dict, exsec_dict, phylostratr_dict, summary_dict):
    species_list = [key for key in exsec_dict.keys()]

    for phylostrata, level in levels_dict.items():
        summary_dict[level] = {species: [] for species in species_list}

    for species, protein_list in exsec_dict.items():
        for protein in protein_list:
            summary_dict[levels_dict[phylostratr_dict[species][protein]]][species].append(protein)


def output_writing(output, exsec_tag, summary_dict, exsec_dict):
    species_list = [key for key in exsec_dict.keys()]

    with open("{output}_potential_{tag}_exsec_proteins_in_phylostrates.tsv".format(
            output=output, tag=exsec_tag), 'a') as output_file:

        output_file.write("Levels\t{species}\n".format(species="\t".join(species_list)))
        for level, species_values in summary_dict.items():
            values = ["{len} ({percent}%)".format(len=len(species_values[species]),
                      percent=round((len(species_values[species])/len(exsec_dict[species]))*100, 2))
                      for species in species_list]
            output_file.write("{level}\t{values}\n".format(level=level, values="\t".join(values)))


if __name__ == "__main__":
    levels_dict, phylostratr_dict, exsec_dict, summary_dict = {}, {}, {}, {}
    levels_parsing(args.levels, levels_dict)
    phylostratr_parsing(args.phylostratr_merged, phylostratr_dict)
    exsec_parsing(args.exsec_merged, exsec_dict)
    generalization(levels_dict, exsec_dict, phylostratr_dict, summary_dict)
    output_writing(args.output, args.exsec_tag, summary_dict, exsec_dict)