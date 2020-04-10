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


def orthogroups_parsing(orthogroups, species_list, tag_dict, ortho_dict):
    #   species_list.extend(orthogroups.readline().strip().split("\t")[:-1])
    header = orthogroups.readline().strip().split("\t")[:-1]
    for species in header:
        if species != "Orthogroup":
            species_list.append(species)
    # print(species_list)
    print("Species: {species_list}".format(species_list="\t".join(species_list)))
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, counts = description[0], description[1:]
        ortho_dict[group_ID] = {tag_dict[species]: int(counts[species_list.index(species)])
                                for species in species_list}
        # print("{group_ID}: {values}".format(group_ID=group_ID, values=ortho_dict[group_ID]))


def summary(tag_dict, ortho_dict, summary_dict):
    for group_ID, counts in ortho_dict.items():
        species_with_nonzero = []
        for species, tag in tag_dict.items():
            if counts[tag] != 0:
                species_with_nonzero.append(tag)
        key = "|".join(species_with_nonzero)
        if key not in summary_dict:
            summary_dict[key] = []

        summary_dict[key].append(group_ID)


def output_writing(output,  summary_dict):
    with open("{output}.common_orthogroups_summary.tsv".format(output=output), 'a') as ortho_summary:
        ortho_summary.write("Species\tNumber_of_orthogroups_with_only_these_species_presented\tCommon_groups_IDs\n")
        for key, groups in summary_dict.items():
            ortho_summary.write("{key}\t{number}\t{groups}\n".format(key=key, number=len(groups), groups=";".join(groups)))


if __name__ == "__main__":
    species_list, tag_dict, ortho_dict, summary_dict = [], {}, {}, {}
    short_name_parsing(args.species_tag, tag_dict)
    orthogroups_parsing(args.ortho, species_list, tag_dict, ortho_dict)
    summary(tag_dict, ortho_dict, summary_dict)
    output_writing(args.output, summary_dict)