try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--summary_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def parsing_summary_table(summary_table, summary_dict):
    header = summary_table.readline()
    for line in summary_table:
        description = line.strip().split("\t")
        gene, trans, prot, ortho, name, nt, nt_identity, nr, nr_identity, swiss, swiss_identity, neuro, \
        neuro_identity, domains, go, ec, ko, pathway, module, reaction, rclass, brite, tc, cazy, bigg, \
        og, cog, eggnog_description, externa, growing, middle, terminal, whole, specificity, markers, exsec = description[0], \
        description[1], description[2], description[3], description[4], description[5], description[6], \
        description[7], description[8], description[9], description[10], description[11], description[12], \
        description[13], description[14], description[15], description[16], description[17], description[18], \
        description[19], description[20], description[21], description[22], description[23], description[24], \
        description[25], description[26], description[27], description[28], description[29], description[30], \
        description[31], description[32], description[33], description[34], description[35]
        if exsec == "Potential_nonclassical" or exsec == "Potential_classical":
            summary_dict[gene] = {"ortho": ortho, "name": name, "nt": nt, "nr": nr, "swiss": swiss,
                                  "neuro": neuro, "domains": domains, "go": go, "ec": ec, "ko": ko,
                                  "pathway": pathway, "module": module, "reaction": reaction,
                                  "rclass": rclass, "brite": brite, "tc": tc, "cazy": cazy,
                                  "bigg": bigg, "og": og, "cog": cog, "eggnog_description": eggnog_description,
                                  "specificity": specificity, "markers": markers, "exsec": exsec}
    print("Length of the summary_dict: {len}".format(len=len(summary_dict.keys())))


def summary(summary_dict, output):
    classical_ortho, nonclassical_ortho = [], []
    classical_specificity, nonclassical_specificity = [], []
    classical_markers, nonclassical_markers = [], []
    classical_undef, nonclassical_undef = [], []

    annotation_keys = ["name", "nt", "nr", "swiss", "neuro", "go", "ec", "ko", "pathway",
                       "module", "reaction", "rclass", "brite", "tc", "cazy", "bigg", "og",
                       "cog", "eggnog_description"]

    for gene, values in summary_dict.items():
        if values["ortho"] != "-":
            if values["exsec"] == "Potential_classical":
                classical_ortho.append(gene)
            elif values["exsec"] == "Potential_nonclassical":
                nonclassical_ortho.append(gene)

        if values["specificity"] != "-":
            if values["exsec"] == "Potential_classical":
                classical_specificity.append(gene)
            elif values["exsec"] == "Potential_nonclassical":
                nonclassical_specificity.append(gene)

        if values["markers"] != "-":
            if values["exsec"] == "Potential_classical":
                classical_markers.append(gene)
            elif values["exsec"] == "Potential_nonclassical":
                nonclassical_markers.append(gene)

        annotation = [values[key] for key in annotation_keys]
        if len(set(annotation)) == 1 and annotation[0] == "-":
            if values["exsec"] == "Potential_classical":
                classical_undef.append(gene)
            elif values["exsec"] == "Potential_nonclassical":
                nonclassical_undef.append(gene)

    with open("{output}.summary.tsv".format(output=output), 'a') as output_file:
        output_file.write("Metrics\tClassical_exsec\tNonclassical_exsec\n")
        output_file.write("Included_in_orthogroups\t{classical}\t{nonclassical}\n".format(
            classical=len(classical_ortho), nonclassical=len(nonclassical_ortho)))
        output_file.write("Preferential_expression\t{classical}\t{nonclassical}\n".format(
            classical=len(classical_specificity), nonclassical=len(nonclassical_specificity)))
        output_file.write("Molecular_markers\t{classical}\t{nonclassical}\n".format(
            classical=len(classical_markers), nonclassical=len(nonclassical_markers)))
        output_file.write("Without_annotations\t{classical}\t{nonclassical}\n".format(
            classical=len(classical_undef), nonclassical=len(nonclassical_undef)))


def output_writing(summary_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tOrthogroup\tEggNOG:Preferred_name\t"
                          "NCBInt\tNCBInr\tSwissProt\tNeuroPep\t"
                          "Domains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tSpecificity\tMolecular_markers\t"
                          "Secretory\Excretory_sequences\n")
        for gene, values in summary_dict.items():
            output_file.write("{gene}\t{ortho}\t{name}\t{nt}\t{nr}\t{swiss}\t{neuro}\t{domains}\t{go}\t{ec}\t"
                              "{ko}\t{pathway}\t{module}\t{reaction}\t{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t"
                              "{og}\t{cog}\t{eggnog_description}\t{specificity}\t{markers}\t{exsec}\n".format(
                                gene=gene, ortho=values["ortho"], name=values["name"], nt=values["nt"],
                                nr=values["nr"], swiss=values["swiss"], neuro=values["neuro"], domains=values["domains"],
                                go=values["go"], ec=values["ec"], ko=values["ko"], pathway=values["pathway"],
                                module=values["module"], reaction=values["reaction"], rclass=values["rclass"],
                                brite=values["brite"], tc=values["tc"], cazy=values["cazy"], bigg=values["bigg"],
                                og=values["og"], cog=values["cog"], eggnog_description=values["eggnog_description"],
                                specificity=values["specificity"], markers=values["markers"], exsec=values["exsec"]
            ))


if __name__ == "__main__":
    summary_dict = {}
    parsing_summary_table(args.summary_table, summary_dict)
    summary(summary_dict, args.output)
    output_writing(summary_dict, args.output)