try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--psi_tab', type=argparse.FileType('r'), required=True,
                    help="Summary table for Psilotrema simillimum")
parser.add_argument('--sps_tab', type=argparse.FileType('r'), required=True,
                    help="Summary table for Sphaeridiotrema pseudoglobulus")
parser.add_argument('--orthogroups_summary', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def summary_table_parsing(species_tab, species_dict):
    header = species_tab.readline()
    for line in species_tab:
        description = line.strip().split("\t")
        protein, orthogroup, preferred_name, nr, swiss, \
        domains, go, ec, ko, pathway, module, reaction, rclass, \
        brite, tc, cazy, bigg, og, cog = description[2], description[3], description[4], description[7], \
                                         description[9], description[11], description[12], description[13], \
                                         description[14], description[15], description[16], description[17], \
                                         description[18], description[19], description[20], description[21], \
                                         description[22], description[23], description[24]
        if orthogroup != "-":
            species_dict[protein] = {"Orthogroup": orthogroup, "EggNOG:Preferred_name": preferred_name,
                                     "NCBInr": nr, "SwissProt": swiss, "Domains_arch": domains,
                                     "EggNOG:GO_terms": go, "EggNOG:EC_number": ec, "EggNOG:KEGG_KO": ko,
                                     "EggNOG:KEGG_Pathway": pathway, "EggNOG:KEGG_Module": module,
                                     "EggNOG:KEGG_Reaction": reaction, "EggNOG:rclass": rclass,
                                     "EggNOG:BRITE": brite, "EggNOG:KEGG_TC": tc, "EggNOG:CAZy": cazy,
                                     "EggNOG:BiGG_Reaction": bigg, "EggNOG:OG": og, "EggNOG:COG_cat": cog}


def orthogroups_parsing(orthogroups_summary, groups_dict):
    header = orthogroups_summary.readline()
    for line in orthogroups_summary:
        description = line.strip().split("\t")
        species, number_of_groups, groups = description[0], description[1], description[2].split(";")
        if species == "Psimillimum":
            groups_dict["Psimillimum"] = groups
        elif species == "Psimillimum|Spseudoglobulus":
            groups_dict["Psimillimum|Spseudoglobulus"] = groups
        elif species == "Spseudoglobulus":
            groups_dict["Spseudoglobulus"] = groups


def append_protein(species, groups, protein, values, output_file):
    if values["Orthogroup"] in groups:
        output_file.write(
            "{ortho}\t{species}\t{protein}\t{name}\t{nr}\t{swiss}\t{domains}\t{go}\t"
            "{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t{rclass}\t{brite}\t{tc}\t"
            "{cazy}\t{bigg}\t{og}\t{cog}\n".format(
                ortho=values["Orthogroup"], species=species, protein=protein,
                name=values["EggNOG:Preferred_name"],
                nr=values["NCBInr"], swiss=values["SwissProt"], domains=values["Domains_arch"],
                go=values["EggNOG:GO_terms"], ec=values["EggNOG:EC_number"],
                ko=values["EggNOG:KEGG_KO"], pathway=values["EggNOG:KEGG_Pathway"],
                module=values["EggNOG:KEGG_Module"], reaction=values["EggNOG:KEGG_Reaction"],
                rclass=values["EggNOG:rclass"], brite=values["EggNOG:BRITE"],
                tc=values["EggNOG:KEGG_TC"], cazy=values["EggNOG:CAZy"],
                bigg=values["EggNOG:BiGG_Reaction"], og=values["EggNOG:OG"],
                cog=values["EggNOG:COG_cat"]))


def output_writing(output, psilo_dict, sphaer_dict, groups_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Orthogroups\tSpecificity\tProteins\tEggNOG:Preferred_names\tNCBInr\tSwissProt\t"
                          "Domains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\tEggNOG:KEGG_KO\t"
                          "EggNOG:KEGG_Pathways\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_categories\n")
        for species, groups in groups_dict.items():
            if species == "Psimillimum":
                for protein, values in psilo_dict.items():
                    append_protein(species, groups, protein, values, output_file)

            elif species == "Spseudoglobulus":
                for protein, values in sphaer_dict.items():
                    append_protein(species, groups, protein, values, output_file)

            elif species == "Psimillimum|Spseudoglobulus":
                for protein, values in psilo_dict.items():
                    append_protein(species, groups, protein, values, output_file)
                for protein, values in sphaer_dict.items():
                    append_protein(species, groups, protein, values, output_file)


if __name__ == "__main__":
    psi_dict, sps_dict, groups_dict = {}, {}, {}
    test_dict = {}
    summary_table_parsing(args.psi_tab, psi_dict)
    summary_table_parsing(args.sps_tab, sps_dict)
    orthogroups_parsing(args.orthogroups_summary, groups_dict)
    output_writing(args.output, psi_dict, sps_dict, groups_dict)

