try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import collections

parser = argparse.ArgumentParser()
parser.add_argument('--orthogroups', type=argparse.FileType('r'), required=True,
                    help="The tsv file with orthogroups (OrthoFinder output file). "
                         "All empty spaces should be replaced with '-'!")
parser.add_argument('--phylostratr_merged', type=argparse.FileType('r'), required=True,
                    help="Merged table with results of phylostratigraphy analysis carried out with Phylostratr."
                         "All species used in orthogroup reconstraction should be in table")
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="The table with all phylostratum (first column) and level description (second column)")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_merged_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        proteinID, mrca, ps, mrca_name = description[0], description[1], description[2], description[3]
        phylostratr_dict[proteinID] = mrca_name


def levels_parsing(levels, levels_dict):
    for line in levels:
        levels_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]


def orthogroups_parsing(orthogroups, ortho_dict, all_species_list):
    all_species_list.extend([species for species in orthogroups.readline().strip().split("\t")[1:]])
    print("Species: {all_species}".format(all_species=",".join(all_species_list)))

    for line in orthogroups:
        description = line.strip().split("\t")
        orthogroup, protein_sets = description[0], description[1:]
        ortho_dict[orthogroup] = {species: protein_sets[all_species_list.index(species)].split(", ")
                                  for species in all_species_list}


def phylostratr_content_of_orthogroups(ortho_dict, phylostratr_dict, levels_dict, all_species_list,
                                       phylostratr_content_dict):
    for species in all_species_list:
        phylostratr_content_dict[species] = {}

    for orthogroup, values in ortho_dict.items():
        for species, protein_set in values.items():
            if len(protein_set) != 1 and protein_set[0] != "-":
                phylostratr_content_dict[species][orthogroup] = {level: [] for level in set(levels_dict.values())}
                for proteinID in protein_set:
                    level_key = levels_dict[phylostratr_dict[proteinID]]
                    phylostratr_content_dict[species][orthogroup][level_key].append(proteinID)


def content_comparison(phylostratr_content_dict, all_species_list, Ratio_dict):
    for species in all_species_list:
        Ratio_dict[species] = {other_species: {} for other_species in all_species_list}

    for first_species in all_species_list:
        for second_species in all_species_list:
            all_common_orthogroups = set.intersection(*[set(phylostratr_content_dict[first_species].keys()),
                                                        set(phylostratr_content_dict[second_species].keys())])
            common_orthogroups_with_same_content = []
            for common_orthogroup in all_common_orthogroups:
                first_species_content = \
                    [level for level in phylostratr_content_dict[first_species][common_orthogroup].keys()
                        if len(phylostratr_content_dict[first_species][common_orthogroup][level]) != 0]

                second_species_content = \
                    [level for level in phylostratr_content_dict[second_species][common_orthogroup].keys()
                     if len(phylostratr_content_dict[second_species][common_orthogroup][level]) != 0]

                if collections.Counter(first_species_content) == collections.Counter(second_species_content):
                    common_orthogroups_with_same_content.append(common_orthogroup)

            Ratio_dict[first_species][second_species] = \
                {"Ratio": round(len(common_orthogroups_with_same_content)/len(all_common_orthogroups), 2),
                 "Same_content": common_orthogroups_with_same_content}


def write_counts_or_ids(levels, species, output, phylostratr_content_dict, tag):
    with open("{output}.{species}.orthogroup_phylostratum_content.protein_{tag}.tsv".format(
            output=output, species=species, tag=tag), 'a') as output_phylostratr_content:
        output_phylostratr_content.write("Orthogroup\t{levels}\n".format(levels="\t".join(levels)))

        for orthogroup, values in phylostratr_content_dict[species].items():
            output_values = []
            if tag == "counts":
                output_values.extend([str(len(values[level])) for level in levels])
            elif tag == "IDs":
                output_values.extend(["{ids}".format(
                    ids="|".join(values[level])) if len(values[level]) != 0 else "-" for level in levels])

            output_phylostratr_content.write("{ortho}\t{output_values}\n".format(
                ortho=orthogroup, output_values="\t".join(output_values)))


def write_ratio(output, all_species_list, Ratio_dict, tag):
    with open("{output}.orthogroups_with_similar_phylostratum_content.{tag}.tsv".format(
            output=output, tag=tag), 'a') as comparison_output:

        comparison_output.write("Species\t{all_species}\n".format(all_species="\t".join(all_species_list)))
        for species in all_species_list:
            output_values = []
            for other_species in all_species_list:
                if tag == "ratio_scores":
                    output_values.append(str(Ratio_dict[species][other_species]["Ratio"]))
                elif tag == "orthogroups_counts":
                    output_values.append(str(len(Ratio_dict[species][other_species]["Same_content"])))

            comparison_output.write("{species}\t{output_values}\n".format(
                species=species, output_values="\t".join(output_values)))


def output_writing(output, all_species_list, levels_dict, phylostratr_content_dict, Ratio_dict):
    levels = [level for level in set(levels_dict.values())]

    for species in all_species_list:
        write_counts_or_ids(levels, species, output, phylostratr_content_dict, "counts")
        write_counts_or_ids(levels, species, output, phylostratr_content_dict, "IDs")

    write_ratio(output, all_species_list, Ratio_dict, "ratio_scores")
    write_ratio(output, all_species_list, Ratio_dict, "orthogroups_counts")

    with open("{output}.orthogroups_with_similar_phylostratum_content.orthogroup_ID.tsv".format(
            output=output), 'a') as IDs_of_common:
        pairs = []
        IDs_of_common.write("Species_pairs\tOrthogroups_IDs\n")
        for species in all_species_list:
            for other_species in all_species_list:
                if species != other_species:
                    pair_key, rev_pair_key = "{first}_{second}".format(first=species, second=other_species), \
                                             "{second}_{first}".format(second=other_species, first=species)
                    if pair_key not in pairs and rev_pair_key not in pairs:
                        IDs_of_common.write("{pair}\t{ids}\n".format(
                            pair=pair_key,
                            ids=";".join(Ratio_dict[species][other_species]["Same_content"])))
                        pairs.extend([pair_key, rev_pair_key])


if __name__ == "__main__":
    all_species_list, ortho_dict, phylostratr_dict, phylostratr_content_dict, levels_dict, Ratio_dict = \
        [], {}, {}, {}, {}, {}
    print("***** Input files parsing *****")
    phylostratr_merged_parsing(args.phylostratr_merged, phylostratr_dict)
    levels_parsing(args.levels, levels_dict)
    orthogroups_parsing(args.orthogroups, ortho_dict, all_species_list)
    print("***** Data analysis *****")
    phylostratr_content_of_orthogroups(ortho_dict, phylostratr_dict, levels_dict, all_species_list,
                                       phylostratr_content_dict)
    content_comparison(phylostratr_content_dict, all_species_list, Ratio_dict)
    print("***** Output files writing *****")
    output_writing(args.output, all_species_list, levels_dict, phylostratr_content_dict, Ratio_dict)