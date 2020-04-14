try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="OrthoFinder output file: Orthogroups.GeneCount.csv")
parser.add_argument('--species_tag', type=argparse.FileType('r'), required=True,
                    help="Table with short name of species (one per line; tab-separated): "
                         "Csinensis_UP000286415.100aa Csinensis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def short_name_parsing(species_tag, tag_dict):
    for line in species_tag:
        description = line.strip().split("\t")
        tag_dict[description[0]] = description[1]


def orthogroups_parsing(orthogroups, species_list, ortho_dict):
    #   species_list.extend(orthogroups.readline().strip().split("\t")[:-1])
    header = orthogroups.readline().strip().split("\t")[:-1]
    for species in header:
        if species != "Orthogroup":
            species_list.append(species)
            ortho_dict[species] = []

    print("Species: {species_list}".format(species_list="\t".join(species_list)))
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, counts = description[0], description[1:]
        for species in ortho_dict.keys():
            if int(counts[species_list.index(species)]) != 0:
                ortho_dict[species].append(group_ID)


def output_writing(output, ortho_dict, tag_dict):
    for species, groups in ortho_dict.items():
        with open("{output}.orthogroups_with_{species}.tsv".format(output=output, species=tag_dict[species]), 'a') \
                as output_for_species:
            output_for_species.write("{groups}".format(groups="\n".join(groups)))


if __name__ == "__main__":
    species_list, ortho_dict, tag_dict = [], {}, {}
    short_name_parsing(args.species_tag, tag_dict)
    orthogroups_parsing(args.ortho, species_list, ortho_dict)
    output_writing(args.output, ortho_dict, tag_dict)