try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="fasta file with proteins")
parser.add_argument('--nr', type=argparse.FileType('r'), required=True)
parser.add_argument('--swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--orthogroups', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def fasta_parsing(protein_dict, fasta):
    protein_fasta = SeqIO.parse(fasta, 'fasta')
    for seq in protein_fasta:
        protein_dict[seq.id] = {"Domains": [], "Anno_swiss": {"hit": [], "identity": []},
                                "Anno_nr": {"hit": [], "identity": []}, "Orthogroup": [],
                                "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                                "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                                "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                                "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                                "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": []}


def orthogroups_parsing(protein_dict, orthogroups):
    head = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(","):
                if protein in protein_dict.keys():
                    protein_dict[protein]["Orthogroup"].append(group_ID)

    for contig, values in protein_dict.items():
        if len(values["Orthogroup"]) == 0:
            values["Orthogroup"].append("-")


def BLAST_annotation(protein_dict, BLAST, key_tag):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        ID, hit_name = description[0].strip(), description[3]

        if ID in protein_dict.keys():
            if hit_name == "no hits":
                protein_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
                protein_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
            else:
                protein_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
                # Per cent identity:
                protein_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append(description[9])


def domains_parsing(protein_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        ID, domain_name = description[0].strip(), description[1]
        protein_dict[ID]["Domains"].append(domain_name)

    for protein in protein_dict.keys():
        if len(protein_dict[protein]["Domains"]) == 0:
            protein_dict[protein]["Domains"].append("-")


def eggNOG_mapper_parsing(protein_dict, eggNOG):
    keys = ["EggNOG:Preferred_name", "EggNOG:GO_terms", "EggNOG:EC_number", "EggNOG:KEGG_KO", "EggNOG:KEGG_Pathway",
            "EggNOG:KEGG_Module", "EggNOG:KEGG_Reaction", "EggNOG:rclass", "EggNOG:BRITE", "EggNOG:KEGG_TC",
            "EggNOG:CAZy", "EggNOG:BiGG_Reaction", "EggNOG:OG", "EggNOG:COG_cat", "EggNOG:Description"]
    index_in_annotation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20]
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein_ID, annotation = description[0], description[1:]
            for el in index_in_annotation:
                if len(annotation) >= el + 1 and len(annotation[el]) != 0:
                    protein_dict[protein_ID][keys[index_in_annotation.index(el)]].append(annotation[el])

    for protein_ID, values in protein_dict.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append('-')


def write_output(protein_dict, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        output.write("### EggNOG:Preferred_name - Predicted gene name\n"
                     "### Orthogroups - ID of the orthogroup constructed for Amoeboaphelidium protococcarum (X4) "
                     "and Paraphelidium tribonemae\n"
                     "### NCBInr_hit - best hit with NCBI non-redundant protein database\n"
                     "### NCBInr_identity(%) - best hit identity percentage\n"
                     "### SwissProt_hit - best hit with SwissProt protein database\n"
                     "### SwissProt_identity - best hit identity percentage\n"
                     "### Domain_arch - reconstructed protein architecture based on analysis of comparison results "
                     "with a PfamAdatabase\n"
                     "### EggNOG:GO_terms - comma delimited list of predicted Gene Ontology terms\n"
                     "### EggNOG:EC_number - the Enzyme Commission number (EC number)\n"
                     "### EggNOG:KEGG_KO - comma delimited list of predicted KEGG KOs\n"
                     "### EggNOG:KEGG_Pathway - comma delimited list of predicted KEGG Pathways\n"
                     "### EggNOG:KEGG_Module - predicted KEGG Module\n"
                     "### EggNOG:KEGG_Reaction - predicted KEGG Reaction\n"
                     "### EggNOG:rclass - predicted KEGG Reaction Class\n"
                     "### EggNOG:BRITE - predicted KEGG BRITE\n"
                     "### EggNOG:KEGG_TC - predicted KEGG TC\n"
                     "### EggNOG:CAZy - best hit with Carbohydrate-Active enZYmes Database\n"
                     "### EggNOG:BiGG_Reaction - comma delimited list of predicted BiGG metabolic reactions\n"
                     "### EggNOG:OG - comma delimited list of matching eggNOG Orthologous Groups\n"
                     "### EggNOG:COG_cat - COG functional category inferred from best matching OG\n"
                     "### EggNOG:Description - eggNOG functional description inferred from best matching OG\n"
                     "# Protein_ID\tEggNOG:Preferred_name\tOrthogroup\tNCBInr_hit\tNCBInr_identity (%)\t"
                     "SwissProt_hit\tSwissProt_identity (%)\tDomain_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                     "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                     "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                     "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\n")
        for protein_ID, values in protein_dict.items():
            output.write("{pep}\t{name}\t{ortho}\t{nr}\t{nr_percent}\t{swiss}\t{swiss_percent}\t"
                         "{domain_arch}\t{go}\t{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t"
                         "{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\n".format(
                          pep=protein_ID, name=values["EggNOG:Preferred_name"][0], ortho=values["Orthogroup"][0],
                          nr=values["Anno_nr"]["hit"][0], nr_percent=values["Anno_nr"]["identity"][0],
                          swiss=values["Anno_swiss"]["hit"][0], swiss_percent=values["Anno_swiss"]["identity"][0],
                          domain_arch="|".join(values["Domains"]), go=values["EggNOG:GO_terms"][0],
                          ec=values["EggNOG:EC_number"][0], ko=values["EggNOG:KEGG_KO"][0],
                          pathway=values["EggNOG:KEGG_Pathway"][0], module=values["EggNOG:KEGG_Module"][0],
                          reaction=values["EggNOG:KEGG_Reaction"][0], rclass=values["EggNOG:rclass"][0],
                          brite=values["EggNOG:BRITE"][0], tc=values["EggNOG:KEGG_TC"][0], cazy=values["EggNOG:CAZy"][0],
                          bigg=values["EggNOG:BiGG_Reaction"][0], og=values["EggNOG:OG"][0],
                          cog=values["EggNOG:COG_cat"][0], description=values["EggNOG:Description"][0]))


if __name__ == "__main__":
    protein_dict = {}
    print("***** Input files parsing *****")
    fasta_parsing(protein_dict, args.fasta)
    orthogroups_parsing(protein_dict, args.orthogroups)
    BLAST_annotation(protein_dict, args.nr, "nr")
    BLAST_annotation(protein_dict, args.swiss, "swiss")
    domains_parsing(protein_dict, args.domains)
    eggNOG_mapper_parsing(protein_dict, args.eggnog)
    print("***** Output file writing *****")
    write_output(protein_dict, args.out)