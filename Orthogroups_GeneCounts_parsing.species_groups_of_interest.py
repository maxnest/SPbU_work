try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="OrthoFinder output file: Orthogroups.GeneCount.csv")
parser.add_argument('--groups', type=argparse.FileType('r'), required=True,
                    help="Table with interesting groups of species (one species per line; tab-separated): "
                         "Csinensis_UP000286415.100aa Redioid_species")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def groups_of_interest_parsing(species_tag, groups_dict):
    for line in species_tag:
        description = line.strip().split("\t")
        if description[1] not in groups_dict.keys():
            groups_dict[description[1]] = []
        groups_dict[description[1]].append(description[0])

    for species_group, species_list in groups_dict.items():
        print("{species_group}: {species_list}".format(species_group=species_group,
                                                       species_list=";".join(species_list)))


def orthogroups_parsing(orthogroups, species_list, ortho_dict):
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
        ortho_dict[group_ID] = {species: int(counts[species_list.index(species)]) for species in species_list}


def summary(group_dict, ortho_dict, summary_dict):
    for species_group, species_list in group_dict.items():
        summary_dict[species_group] = {}
        if species_group == "Psimillimum" or species_group == "Spseudoglobulus":
            for group_ID, values in ortho_dict.items():
                for species in species_list:
                    if values[species] != 0:
                        values_in_redioid = [values[species] for species in groups_dict["Redioid_species"]]
                        values_in_sporocystoid = [values[species] for species in groups_dict["Sporocystoid_species"]]
                        values_in_freeliving = [values[species] for species in groups_dict["Freeliving_species"]]
                        if values_in_redioid.count(0) == 0 \
                                or values_in_sporocystoid.count(0) == 0 \
                                or values_in_freeliving.count(0) == 0:
                            summary_dict[species_group][group_ID] = {species: values[species]
                                                                     for species in species_list}

                        elif values_in_redioid.count(0) == len(values_in_redioid) \
                                and values_in_sporocystoid.count(0) == len(values_in_sporocystoid) \
                                and values_in_freeliving.count(0) == len(values_in_freeliving):
                            summary_dict[species_group][group_ID] = {species: values[species]
                                                                     for species in species_list}
        else:
            for group_ID, values in ortho_dict.items():
                species_with_nonzero = []
                for species in species_list:
                    if values[species] != 0:
                        species_with_nonzero.append(species)
                if len(species_with_nonzero) == len(species_list):
                    summary_dict[species_group][group_ID] = {species: values[species] for species in species_list}


def output_writing(output,  summary_dict, groups_dict):
    for species_group in summary_dict.keys():
        with open("{output}.{species_group}.tsv".format(output=output,
                                                        species_group=species_group), 'a') as output_file:
            output_file.write("Group_ID\t{species}\n".format(species="\t".join(groups_dict[species_group])))
            for group_ID, values in summary_dict[species_group].items():
                output_file.write("{group_ID}\t{values}\n".format(group_ID=group_ID,
                                                                  values="\t".join([str(values[species]) for species in
                                                                                    groups_dict[species_group]])))


if __name__ == "__main__":
    species_list, groups_dict, ortho_dict, summary_dict = [], {}, {}, {}
    groups_of_interest_parsing(args.groups, groups_dict)
    orthogroups_parsing(args.ortho, species_list, ortho_dict)
    summary(groups_dict, ortho_dict, summary_dict)
    output_writing(args.output, summary_dict, groups_dict)
