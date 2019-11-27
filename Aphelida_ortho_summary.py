try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="Orthogroups.csv")
parser.add_argument('--unassigned', type=argparse.FileType('r'), required=True,
                    help="Orthogroups_UnassignedGenes.csv")
parser.add_argument('--x4_anno', type=argparse.FileType('r'), required=True,
                    help="eggNOG-mapper output file for X4")
parser.add_argument('--para_anno', type=argparse.FileType('r'), required=True,
                    help="eggNOG-mapper output file for Paraphelidium")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def orthogroup_parsing(ortho, ortho_dict):
    head = ortho.readline()
    for line in ortho:
        description = line.strip().split("\t")
        group_ID, proteins_list = description[0], description[1:]
        ortho_dict[group_ID] = {"X4": [], "Para": [], "COGs": [], "Description": []}
        for proteins in proteins_list:
            for protein in proteins.split(","):
                if "Aprotococcarum_X4" in protein:
                    ortho_dict[group_ID]["X4"].append(protein)
                elif "Partr_" in protein:
                    ortho_dict[group_ID]["Para"].append(protein)

        if len(ortho_dict[group_ID]["X4"]) != 0 and len(ortho_dict[group_ID]["Para"]) != 0:
            ortho_dict[group_ID]["Description"].append("Common")
        elif len(ortho_dict[group_ID]["X4"]) != 0 and len(ortho_dict[group_ID]["Para"]) == 0:
            ortho_dict[group_ID]["Description"].append("X4-specific")
        elif len(ortho_dict[group_ID]["X4"]) == 0 and len(ortho_dict[group_ID]["Para"]) != 0:
            ortho_dict[group_ID]["Description"].append("Para-specific")


def unassigned_parsing(unassigned, x4_unassigned, para_unassigned):
    head = unassigned.readline()
    for line in unassigned:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein in proteins:
            if "Aprotococcarum_X4" in protein:
                x4_unassigned.append(protein)
            elif "Partr_" in protein:
                para_unassigned.append(protein)


def anno_parsing(anno, anno_dict):
    for line in anno:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            pep_ID, annotation = description[0], description[1:]
            if len(annotation) >= 20:
                COGs = list(annotation[19])  # one gene may have several annotations
                anno_dict[pep_ID] = COGs


def orthogroup_anno_analysis(ortho_dict, x4_anno_dict, para_anno_dict):
    for group, values in ortho_dict.items():
        annotations = []
        if values["Description"][0] == "Common":
            for protein in values["X4"]:
                if protein in x4_anno_dict.keys():
                    annotations.extend(x4_anno_dict[protein])

            for protein in values["Para"]:
                if protein in para_anno_dict.keys():
                    annotations.extend(para_anno_dict[protein])

        elif values["Description"][0] == "X4-specific":
            for protein in values["X4"]:
                if protein in x4_anno_dict.keys():
                    annotations.extend(x4_anno_dict[protein])

        elif values["Description"][0] == "Para-specific":
            for protein in values["Para"]:
                if protein in para_anno_dict.keys():
                    annotations.extend(para_anno_dict[protein])

        if len(annotations) != 0:
            # we only exclude repetitions, and do not search for the most frequent,
            # since there are few or no cases when there is only 1 annotation with a maximum frequency
            values["COGs"].extend(set(annotations))


def summary(out, ortho_dict, x4_unassigned, para_unassigned, x4_anno_dict, para_anno_dict):
    COG_dict = {}
    common_in_ortho, x4_specific_in_ortho, para_specific_in_ortho = [], [], []
    common_in_COG, x4_specific_in_COG, para_specific_in_COG, total_with_COGs = [], [], [], []
    COG_description = ["D:Cell cycle control, cell division, chromosome partitioning",
                       "M:Cell wall/membrane/envelope biogenesis", "N:Cell motility",
                       "O:Post-translational modification, protein turnover, and chaperones",
                       "T:Signal transduction mechanisms",
                       "U:Intracellular trafficking, secretion, and vesicular transport",
                       "V:Defense mechanisms", "W:Extracellular structures",
                       "Y:Nuclear structure", "Z:Cytoskeleton", "A:RNA processing and modification",
                       "B:Chromatin structure and dynamics", "J:Translation, ribosomal structure and biogenesis",
                       "K:Transcription", "L:Replication, recombination and repair",
                       "C:Energy production and conversion", "E:Amino acid transport and metabolism",
                       "F:Nucleotide transport and metabolism", "G:Carbohydrate transport and metabolism",
                       "H:Coenzyme transport and metabolism", "I:Lipid transport and metabolism",
                       "P:Inorganic ion transport and metabolism",
                       "Q:Secondary metabolites biosynthesis, transport, and catabolism",
                       "R:General function prediction only", "S:Function unknown"]
    for COG in COG_description:
        COG_dict[COG.split(":")[0]] = {"Description": COG.split(":")[1], "Common_groups": [], "X4_specific_groups": [],
                                       "Para_specific_groups": [], "X4_unassigned": [], "Para_unassigned": []}

    for group, values in ortho_dict.items():
        if values["Description"][0] == "Common":
            common_in_ortho.append(group)
        elif values["Description"][0] == "X4-specific":
            x4_specific_in_ortho.append(group)
        elif values["Description"][0] == "Para-specific":
            para_specific_in_ortho.append(group)

        if len(values["COGs"]) != 0:
            for COG in values["COGs"]:
                if values["Description"][0] == "Common":
                    COG_dict[COG]["Common_groups"].append(group)
                    common_in_COG.append(group)
                    total_with_COGs.append(group)
                elif values["Description"][0] == "X4-specific":
                    COG_dict[COG]["X4_specific_groups"].append(group)
                    x4_specific_in_COG.append(group)
                    total_with_COGs.append(group)
                elif values["Description"][0] == "Para-specific":
                    COG_dict[COG]["Para_specific_groups"].append(group)
                    para_specific_in_COG.append(group)
                    total_with_COGs.append(group)

    for protein in x4_unassigned:
        if protein in x4_anno_dict.keys():
            for COG in x4_anno_dict[protein]:
                COG_dict[COG]["X4_unassigned"].append(protein)

    for protein in para_unassigned:
        if protein in para_anno_dict.keys():
            for COG in para_anno_dict[protein]:
                COG_dict[COG]["Para_unassigned"].append(protein)

    with open("{out}.COGs_summary.tsv".format(out=out), 'a') as output:
        output.write("### Number of orthogroups (# with annotation): {orthogroups_counts} ({orthogroups_anno})\n"
                     "### Number of groups with both species (# with annotation): {both} ({both_anno})\n"
                     "### Number of A.protococcarum (X4) specific groups (# with annotation): {x4_groups} ({x4_anno})\n"
                     "### Number of P.tribonemae specific groups (# with annotation): {para_groups} ({para_anno})\n"
                     "### Number of A.protococcarum (X4) unassigned proteins: {x4_unassigned}\n"
                     "### Number of P.tribonemae unassigned proteins: {para_unassigned}\n"
                     "# 'COG_in_common_groups' - number of common for both species groups where this COG is present\n"
                     "# 'COG_in_X4_specific' - number of A.protococcarum-specific groups where this COG is present\n"
                     "# 'COG_in_Para_specific' - number of P.tribonemae-specific groups where this COG is present\n"
                     "# 'Counts_in_X4_unassigned' - number of A.protococcarum unassigned proteins with this COG\n"
                     "# 'Counts_in_Para_unassigned' - number of P.tribonemae unassigned proteins with this COG\n"
                     "COG\tDescription\tCOG_in_common_groups\tCOG_in_X4_specific\tCOG_in_Para_specific\t"
                     "Counts_in_X4_unassigned\tCounts_in_Para_unassigned\n".format(
                      orthogroups_counts=len(ortho_dict.keys()), orthogroups_anno=len(set(total_with_COGs)),
                      both=len(common_in_ortho), both_anno=len(set(common_in_COG)),
                      x4_groups=len(x4_specific_in_ortho), x4_anno=len(set(x4_specific_in_COG)),
                      para_groups=len(para_specific_in_ortho), para_anno=len(set(para_specific_in_COG)),
                      x4_unassigned=len(x4_unassigned), para_unassigned=len(para_unassigned)))
        for group, values in COG_dict.items():
            output.write("{COG}\t{Description}\t{Common}\t{X4_specific}\t{Para_specific}\t"
                         "{X4_unassigned}\t{Para_unassigned}\n".format(
                          COG=group, Description=values["Description"], Common=len(values["Common_groups"]),
                          X4_specific=len(values["X4_specific_groups"]),
                          Para_specific=len(values["Para_specific_groups"]),
                          X4_unassigned=len(values["X4_unassigned"]), Para_unassigned=len(values["Para_unassigned"])))


if __name__ == "__main__":
    ortho_dict, x4_unassigned, para_unassigned = {}, [], []
    x4_anno_dict, para_anno_dict = {}, {}
    print("***** OrthoFinder output files parsing *****")
    orthogroup_parsing(args.ortho, ortho_dict)
    unassigned_parsing(args.unassigned, x4_unassigned, para_unassigned)
    print("***** eggNOG-mapper output files parsing *****")
    anno_parsing(args.x4_anno, x4_anno_dict)
    anno_parsing(args.para_anno, para_anno_dict)
    print("***** Analysis of orthogroups *****")
    orthogroup_anno_analysis(ortho_dict, x4_anno_dict, para_anno_dict)
    print("***** Summarizing and output file creating *****")
    summary(args.out, ortho_dict, x4_unassigned, para_unassigned, x4_anno_dict, para_anno_dict)


