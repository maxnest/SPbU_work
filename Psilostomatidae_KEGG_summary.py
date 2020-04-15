try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--summary_table', type=argparse.FileType('r'), required=True,
                    help="Summary table with results of analyses: good.*_ref.genes_level.summary_table.tsv")
parser.add_argument('--kegg_description', type=argparse.FileType('r'), required=True,
                    help="Table with description of KEGG pathways: Super-pathway(Category), Code, Description")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def kegg_parsing(kegg_description, kegg_dict):
    header = kegg_description.readline()
    for line in kegg_description:
        description = line.strip().split("\t")
        category, code, pathway_description = description[0], description[1], description[2]
        if category not in kegg_dict.keys():
            kegg_dict[category] = {}

        kegg_dict[category][code] = {"Description": pathway_description,
                                     "Redia": [], "Cercaria": [], "Marita": [], "Total_count": []}


def summary_table_parsing(summary_table, kegg_dict):
    header = summary_table.readline()
    for line in summary_table:
        description = line.strip().split("\t")
        gene, kegg_pathways, redia_exp, cercaria_exp, marita_exp = description[0], description[15], \
            float(description[26]), float(description[27]), float(description[28])

        for kegg_pathway in kegg_pathways.split(","):
            if kegg_pathway.startswith("ko"):
                for category, codes in kegg_dict.items():
                    if kegg_pathway.split("ko")[1] in codes:
                        codes[kegg_pathway.split("ko")[1]]["Total_count"].append(gene)

                        if redia_exp >= 1:
                            codes[kegg_pathway.split("ko")[1]]["Redia"].append(gene)

                        if cercaria_exp >= 1:
                            codes[kegg_pathway.split("ko")[1]]["Cercaria"].append(gene)

                        if marita_exp >= 1:
                            codes[kegg_pathway.split("ko")[1]]["Marita"].append(gene)


def output_writing(output, kegg_dict):
    with open("{output}.all_pathways_with_5_or_more_genes.tsv".format(output=output), 'a') as all_active_pathways:
        all_active_pathways.write("Category\tPathway\tRedia:GeneCount\tCercaria:GeneCount\tMarita:GeneCount\t"
                                  "Total_GeneCount\tRedia:Percent\tCercaria:Percent\tMarita:Percent\n")
        for category, codes in kegg_dict.items():
            for code, values in codes.items():
                if len(values["Total_count"]) >= 5:
                    all_active_pathways.write("{category}\t{pathway}\t{redia_count}\t{cercaria_count}\t{marita_count}\t"
                                              "{total}\t{redia_percent}\t{cercaria_percent}\t{marita_percent}\n".format(
                                                category=category,
                                                pathway="{code}|{description}".format(code=code,
                                                                                      description=values["Description"]),
                                                redia_count=len(values["Redia"]),
                                                cercaria_count=len(values["Cercaria"]),
                                                marita_count=len(values["Marita"]),
                                                total=len(values["Total_count"]),
                                                redia_percent=
                                                round((len(values["Redia"])/len(values["Total_count"])) * 100, 2),
                                                cercaria_percent=
                                                round((len(values["Cercaria"])/len(values["Total_count"])) * 100, 2),
                                                marita_percent=
                                                round((len(values["Marita"])/len(values["Total_count"])) * 100, 2)))


if __name__ == "__main__":
    kegg_dict = {}
    kegg_parsing(args.kegg_description, kegg_dict)
    summary_table_parsing(args.summary_table, kegg_dict)
    output_writing(args.output, kegg_dict)
