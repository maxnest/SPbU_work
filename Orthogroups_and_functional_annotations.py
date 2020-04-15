try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--ortho_summary', type=argparse.FileType('r'), required=True,
                    help="Summary table with description species in orthogroups obtained with "
                         "Orthogroups_GeneCounts_parsing.py")
parser.add_argument('--species_groups', type=argparse.FileType('r'), required=True,
                    help="Table with 2 columns: "
                         "1) Group of interest (for instance: Psimillimum|Spseudoglobulus), "
                         "2) Tag with description for group (for instance: Psilostomatidae). One groups per line")
parser.add_argument('--goterms', type=argparse.FileType('r'), required=True,
                    help="Table with GOterms for orthogroups created with KinFin")
parser.add_argument('--pfama', type=argparse.FileType('r'), required=True,
                    help="Table with domains descriptions for orthogroups created with KinFin")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def ortho_summary_parsing(ortho_summary, ortho_dict):
    header = ortho_summary.readline()
    for line in ortho_summary:
        description = line.strip().split("\t")
        species, groups = description[0], description[2]
        ortho_dict[species] = groups.split(";")


def species_groups_parsing(species_groups, groups_of_interest):
    for line in species_groups:
        description = line.strip().split("\t")
        group, tag = description[0], description[1]
        groups_of_interest[group] = tag


def functional_annotation_parsing(annotation, annotation_dict, annotation_tag):
    header = annotation.readline()
    for line in annotation:
        description = line.strip().split("\t")
        orthogroup, domain_ids, domain_description = description[0], description[3], description[4]
        if orthogroup not in annotation_dict.keys():
            annotation_dict[orthogroup] = {"GO": {}, "PfamA": {}}

        annotation_dict[orthogroup][annotation_tag]["Domain_IDs"] = domain_ids
        annotation_dict[orthogroup][annotation_tag]["Domain_description"] = domain_description


def unite_results(ortho_dict, group_of_interest, annotation_dict, united_results):
    for species_group, tag in group_of_interest.items():
        for orthogroup in ortho_dict[species_group]:
            united_results[orthogroup] = {
                "Tag": tag, "GO:Domain_IDs": annotation_dict[orthogroup]["GO"]["Domain_IDs"],
                "GO:Domain_description": annotation_dict[orthogroup]["GO"]["Domain_description"],
                "PfamA:Domain_IDs": annotation_dict[orthogroup]["PfamA"]["Domain_IDs"],
                "PfamA:Domain_description": annotation_dict[orthogroup]["PfamA"]["Domain_description"]}


def output_writing(output, united_results):
    with open("{output}.tsv".format(output=output), 'a') as  output_file:
        output_file.write("Orthogroup_ID\tClassification\tGeneOntology:Domain_IDs\tGeneOntology:Domain_description\t"
                          "PfamA:Domain_IDs\tPfamA:Domain_description\n")
        for orthogroup, values in united_results.items():
            output_file.write("{ortho}\t{classification}\t{go_domains}\t{go_description}\t"
                              "{pfama_domains}\t{pfama_description}\n".format(
                               ortho=orthogroup, classification=values["Tag"], go_domains=values["GO:Domain_IDs"],
                               go_description=values["GO:Domain_description"], pfama_domains=values["PfamA:Domain_IDs"],
                               pfama_description=values["PfamA:Domain_description"]))


if __name__ == "__main__":
    ortho_dict, groups_of_interest, annotation_dict, united_results = {}, {}, {}, {}
    ortho_summary_parsing(args.ortho_summary, ortho_dict)
    species_groups_parsing(args.species_groups, groups_of_interest)
    functional_annotation_parsing(args.goterms, annotation_dict, "GO")
    functional_annotation_parsing(args.pfama, annotation_dict, "PfamA")
    unite_results(ortho_dict, groups_of_interest, annotation_dict, united_results)
    output_writing(args.output, united_results)
