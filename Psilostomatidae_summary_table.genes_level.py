try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--gene_map', type=argparse.FileType('r'), required=True,
                    help="Table with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nt', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nr', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--scaled_tpm', type=argparse.FileType('r'), required=True)
parser.add_argument('--rediae_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--cercariae_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--marita_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def gene_map_parsing(gene_dict, gene_map):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_ID, transcript_ID, protein_ID = description[0], description[1], description[2]
        gene_dict[gene_ID] = {"transcript": transcript_ID, "protein": protein_ID,
                              "domains": [], "Anno_swiss": {"hit": [], "identity": []},
                              "Anno_nt": {"hit": [], "identity": []},
                              "Anno_nr": {"hit": [], "identity": []},
                              "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                              "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                              "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                              "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                              "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                              "Rediae_scaledTPM": 0, "Cercariae_scaledTPM": 0, "Marita_scaledTPM": 0,
                              "Specificity": [], "Cluster": [], "Ortho": []}


def orthogroups_parsing(gene_dict, prot_2_gene_dict, orthogroups):
    header = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(", "):
                if protein in prot_2_gene_dict.keys():
                    gene_dict[prot_2_gene_dict[protein]]["Ortho"].append(group_ID)

    for gene, values in gene_dict.items():
        if len(values["Ortho"]) == 0:
            values["Ortho"].append("-")


def eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, eggNOG):
    keys = ["EggNOG:Preferred_name", "EggNOG:GO_terms", "EggNOG:EC_number", "EggNOG:KEGG_KO", "EggNOG:KEGG_Pathway",
            "EggNOG:KEGG_Module", "EggNOG:KEGG_Reaction", "EggNOG:rclass", "EggNOG:BRITE", "EggNOG:KEGG_TC",
            "EggNOG:CAZy", "EggNOG:BiGG_Reaction", "EggNOG:OG", "EggNOG:COG_cat", "EggNOG:Description"]
    index_in_annotation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20]

    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein, annotation = description[0], description[1:]
            for el in index_in_annotation:
                if len(annotation) >= el + 1 and len(annotation[el]) != 0:
                    if protein in prot_2_gene_dict.keys():
                        gene_dict[prot_2_gene_dict[protein]][keys[index_in_annotation.index(el)]].append(annotation[el])

    for gene, values in gene_dict.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append('-')


def add_match(description, ID, hit_name, gene_dict, seq_dict, key_tag):
    if hit_name == "no hits":
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
    else:
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["identity"].append(description[9])


def BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, BLAST, key_tag):
    header = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        ID, hit_name = description[0], description[3]
        if ID in trans_2_gene_dict.keys():
            add_match(description, ID, hit_name, gene_dict, trans_2_gene_dict, key_tag)
        elif ID in prot_2_gene_dict.keys():
            add_match(description, ID, hit_name, gene_dict, prot_2_gene_dict, key_tag)


def domains_parsing(gene_dict, prot_2_gene_dict, domains):
    header = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        protein, domain_name = description[0].strip(), description[1]
        if protein in prot_2_gene_dict.keys():
            gene_dict[prot_2_gene_dict[protein]]["domains"].append(domain_name)

    for gene in gene_dict.keys():
        if len(gene_dict[gene]["domains"]) == 0:
            gene_dict[gene]["domains"].append("-")


def clusters(gene_dict, table, dict):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)

    for gene in gene_dict.keys():
        for cluster, genes in dict.items():
            if gene in genes:
                gene_dict[gene]["Cluster"].append(cluster)

    for gene, values in gene_dict.items():
        if len(values["Cluster"]) == 0:
            values["Cluster"].append("-")


def specificity(gene_dict, table, tag):
    for line in table:
        gene = line.strip()
        gene_dict[gene]["Specificity"].append(tag)


def expression(gene_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, red, cer, mar = description[0], float(description[1]), float(description[2]), float(description[3])
        gene_dict[gene]["Rediae_scaledTPM"] += red
        gene_dict[gene]["Cercariae_scaledTPM"] += cer
        gene_dict[gene]["Marita_scaledTPM"] += mar


def write_output_files(gene_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tTranscript_ID\tProtein_ID\tOrthogroup\tEggNOG:Preferred_name\t"
                          "NCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\tDomains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tRediae_scaledTPM\tCercariae_scaledTPM\t"
                          "Marita_scaledTPM\tSpecificity\tCluster\n")
        for gene, values in gene_dict.items():
            output_file.write("{gene}\t{transcript}\t{protein}\t{ortho}\t{name}\t{nt}\t{nt_identity}\t{nr}\t"
                              "{nr_identity}\t{swiss}\t"
                              "{swiss_identity}\t{domains}\t{go}\t{ec}\t{ko}\t{kegg}\t{module}\t{reaction}\t"
                              "{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\t"
                              "{red_tpm}\t{cer_tpm}\t{mar_tpm}\t{specificity}\t"
                              "{cluster}\n".format(
                                gene=gene, transcript=values["transcript"],
                                protein=values["protein"], ortho=values["Ortho"][0],
                                name=values["EggNOG:Preferred_name"][0],
                                nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                                nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                                swiss=values["Anno_swiss"]["hit"][0],
                                swiss_identity=values["Anno_swiss"]["identity"][0],
                                domains="|".join(values["domains"]),
                                go=values["EggNOG:GO_terms"][0], ec=values["EggNOG:EC_number"][0],
                                ko=values["EggNOG:KEGG_KO"][0], kegg=values["EggNOG:KEGG_Pathway"][0],
                                module=values["EggNOG:KEGG_Module"][0], reaction=values["EggNOG:KEGG_Reaction"][0],
                                rclass=values["EggNOG:rclass"][0], brite=values["EggNOG:BRITE"][0],
                                tc=values["EggNOG:KEGG_TC"][0], cazy=values["EggNOG:CAZy"][0],
                                bigg=values["EggNOG:BiGG_Reaction"][0], og=values["EggNOG:OG"][0],
                                cog=values["EggNOG:COG_cat"][0], description=values["EggNOG:Description"][0],
                                red_tpm=values["Rediae_scaledTPM"], cer_tpm=values["Cercariae_scaledTPM"],
                                mar_tpm=values["Marita_scaledTPM"], specificity=values["Specificity"][0],
                                cluster=values["Cluster"][0]))


if __name__ == "__main__":
    gene_dict, cluster_dict = {}, {}
    gene_map_parsing(gene_dict, args.gene_map)
    prot_2_gene_dict = {gene_dict[gene]["protein"]: gene for gene in gene_dict.keys()}
    trans_2_gene_dict = {gene_dict[gene]["transcript"]: gene for gene in gene_dict.keys()}
    print("***** Orthogroups parsing *****")
    orthogroups_parsing(gene_dict, prot_2_gene_dict, args.ortho)
    print("***** eggNOG-mapper output parsing *****")
    eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, args.eggnog)
    print("***** BLAST results parsing *****")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_swiss, "swiss")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nt, "nt")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nr, "nr")
    print("***** Domain architecture parsing *****")
    domains_parsing(gene_dict, prot_2_gene_dict, args.domains)
    print("***** Cluster parsing *****")
    clusters(gene_dict, args.clusters, cluster_dict)
    print("***** Files with IDs of stage-specific sequences parsing *****")
    specificity(gene_dict, args.rediae_specific, "R")
    specificity(gene_dict, args.cercariae_specific, "C")
    specificity(gene_dict, args.marita_specific, "M")
    for gene, values in gene_dict.items():
        if len(values["Specificity"]) == 0:
            values["Specificity"].append('H')
    print("***** Table with averaged expression values parsing *****")
    expression(gene_dict, args.scaled_tpm)
    print("***** Output file creating *****")
    write_output_files(gene_dict, args.output)