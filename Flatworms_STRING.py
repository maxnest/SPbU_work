try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphy analysis for species")
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="Table with description of phylostratigraphic levels")
parser.add_argument('--rbbh', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits (RBBH) for species versus STRING database")
parser.add_argument('--string_links', type=argparse.FileType('r'), required=True,
                    help="STRING table with proteins links description")
parser.add_argument('--string_actions', type=argparse.FileType('r'), required=True,
                    help="STRING table with protein actions description")
parser.add_argument('--flatworm_tag', type=str, required=True)
parser.add_argument('--string_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0], description[1], description[2], description[3]
        phylostratr_dict[protein_ID[1:-1]] = mrca_name[1:-1]


def levels_parsing(levels, levels_dict):
    for line in levels:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def rbbh_parsing(rbbh, rbbh_dict,):
    header = rbbh.readline()
    for line in rbbh:
        description = line.strip().split("\t")
        rbbh_pair_ID, flatworm_seq_ID, string_seq_ID = description[0], description[1], description[2]
        rbbh_dict[rbbh_pair_ID] = {"flatworm": flatworm_seq_ID, "string": string_seq_ID}


def string_links_parsing(string_links, string_links_dict):
    header = string_links.readline()
    for line in string_links:
        description = line.strip().split(" ")
        protein_1, protein_2 = description[0], description[1]
        string_links_dict["{protein_1}_{protein_2}".format(protein_1=protein_1, protein_2=protein_2)] = {
            "first_node": protein_1, "second_node": protein_2,
            "flatworms_2_first_node": [], "flatworms_2_second_node": [], "action_modes": [],
            "first_rbbh_pair": [], "second_rbbh_pair": []
        }


def string_actions_parsing(string_actions, string_links_dict):
    header = string_actions.readline()
    for line in string_actions:
        description = line.strip().split("\t")
        protein_1, protein_2, action_mode, is_directional, a_is_acting, score = description[0], description[1], \
                                                                                description[2], description[4], \
                                                                                description[5], description[6]
        pair_key = "{protein_1}_{protein_2}".format(protein_1=protein_1, protein_2=protein_2)
        string_links_dict[pair_key]["action_modes"].append("{mode} (is_directional:{is_directional};"
                                                           "first_protein_is_acting:{a_is_acting})".format(
                                                            mode=action_mode, is_directional=is_directional,
                                                            a_is_acting=a_is_acting))

    for pair, values in string_links_dict.items():
        if len(values["action_modes"]) == 0:
                values["action_modes"].append("-")


def append_flatworms_proteins(string_links_dict, rbbh_dict):
    for pair, values in string_links_dict.items():
        for rbbh_pair, rbbh_values in rbbh_dict.items():
            if rbbh_values["string"] == values["first_node"] and \
                    rbbh_values["flatworm"] not in values["flatworms_2_first_node"]:
                values["flatworms_2_first_node"].append(rbbh_values["flatworm"])
                values["first_rbbh_pair"].append(rbbh_pair)

            if rbbh_values["string"] == values["second_node"] and \
                    rbbh_values["flatworm"] not in values["flatworms_2_second_node"]:
                values["flatworms_2_second_node"].append(rbbh_values["flatworm"])
                values["second_rbbh_pair"].append(rbbh_pair)


def output_writing(output, phylostratr_dict, levels_dict, string_links_dict, flatworm_tag, string_tag):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("{flatworm}_protein_1_ID\tFirst_RBBH_ID\t{flatworm}_protein_1_phylostrata\t"
                          "{flatworm}_protein_2_ID\tSecond_RBBH_ID\t{flatworm}_protein_2_phylostrata\t"
                          "{string}_node_1\tInteraction_type\t{string}_node_2\t"
                          "Interaction_modes_and_directions\n".format(
                            flatworm=flatworm_tag, string=string_tag))

        for pair, values in string_links_dict.items():
            if len(values["flatworms_2_first_node"]) != 0 and len(values["flatworms_2_second_node"]) != 0:
                for first_protein in set(values["flatworms_2_first_node"]):
                    for second_protein in set(values["flatworms_2_second_node"]):
                        output_file.write("{first_protein}\t{first_rbbh}\t{first_phylostrata}\t"
                                          "{second_protein}\t{second_rbbh}\t{second_phylostrata}\t"
                                          "{first_node}\tprotein_2_protein\t{second_node}\t{actions}\n".format(
                                            first_protein=first_protein,
                                            first_rbbh=values["first_rbbh_pair"][
                                                values["flatworms_2_first_node"].index(first_protein)],
                                            first_phylostrata=levels_dict[phylostratr_dict[first_protein]],
                                            second_protein=second_protein,
                                            second_rbbh=values["second_rbbh_pair"][
                                                values["flatworms_2_second_node"].index(second_protein)],
                                            second_phylostrata=levels_dict[phylostratr_dict[second_protein]],
                                            first_node=values["first_node"], second_node=values["second_node"],
                                            actions="|".join(values["action_modes"])))


if __name__ == "__main__":
    phylostratr_dict, levels_dict, rbbh_dict, string_links_dict = {}, {}, {}, {}
    print("***** Input files parsing *****")
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    levels_parsing(args.levels, levels_dict)
    rbbh_parsing(args.rbbh, rbbh_dict)
    string_links_parsing(args.string_links, string_links_dict)
    string_actions_parsing(args.string_actions, string_links_dict)
    print("***** Flatworms proteins appending *****")
    append_flatworms_proteins(string_links_dict, rbbh_dict)
    print("***** Output file writing *****")
    output_writing(args.output, phylostratr_dict, levels_dict, string_links_dict, args.flatworm_tag, args.string_tag)