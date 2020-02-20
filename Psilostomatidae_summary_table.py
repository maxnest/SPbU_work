try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
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


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ID": [seq.id.split(" ")[0]], "nucl": seq.seq, "amino": "",
                                             "domains": [], "Anno_swiss": {"hit": [], "identity": []},
                                             "Anno_nt": {"hit": [], "identity": []},
                                             "Anno_nr": {"hit": [], "identity": []},
                                             "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                                             "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                                             "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                                             "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                                             "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                                             "Rediae_scaledTPM": 0, "Cercariae_scaledTPM": 0, "Marita_scaledTPM": 0,
                                             "Rediae_specific": [], "Cercaria_specific": [], "Marita_specific": [],
                                             "Cluster": [], "Ortho": []}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".p")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"


def orthogroups_parsing(contig_dict, orthogroups):
    head = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(","):
                if protein.split(".p")[0] in contig_dict.keys():
                    contig_dict[protein.split(".p")[0]]["Ortho"].append(group_ID)

    for contig, values in contig_dict.items():
        if len(values["Ortho"]) == 0:
            values["Ortho"].append("-")


def eggNOG_mapper_parsing(contig_dict, eggNOG):
    keys = ["EggNOG:Preferred_name", "EggNOG:GO_terms", "EggNOG:EC_number", "EggNOG:KEGG_KO", "EggNOG:KEGG_Pathway",
            "EggNOG:KEGG_Module", "EggNOG:KEGG_Reaction", "EggNOG:rclass", "EggNOG:BRITE", "EggNOG:KEGG_TC",
            "EggNOG:CAZy", "EggNOG:BiGG_Reaction", "EggNOG:OG", "EggNOG:COG_cat", "EggNOG:Description"]
    index_in_annotation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20]
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            contig_ID, annotation = description[0].split(".p")[0], description[1:]
            for el in index_in_annotation:
                if len(annotation) >= el + 1 and len(annotation[el]) != 0:
                    contig_dict[contig_ID][keys[index_in_annotation.index(el)]].append(annotation[el])

    for contig_ID, values in contig_dict.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append('-')


def BLAST_annotation(contig_dict, BLAST, key_tag):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        if key_tag == "swiss" or key_tag == "nr":
            ID, hit_name = description[0].strip().split(".p")[0], description[3]
        else:
            ID, hit_name = description[0].strip().split(" ")[0], description[3]

        if ID in contig_dict.keys():
            if hit_name == "no hits":
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
            else:
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append(description[9])

    if key_tag == "swiss" or key_tag == "nr":
        for contig in contig_dict.keys():
            if contig_dict[contig]["amino"] == "without ORF":
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["hit"].append("-")
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["identity"].append("-")


def domains_parsing(contig_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        ID, domain_name = description[0].strip().split(".p")[0], description[1]
        contig_dict[ID]["domains"].append(domain_name)

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["domains"]) == 0:
            contig_dict[contig]["domains"].append("-")


def clusters(contig_dict, table, dict):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)

    for contig in contig_dict.keys():
        for cluster, genes in dict.items():
            if contig in genes:
                contig_dict[contig]["Cluster"].append(cluster)

    for contig, values in contig_dict.items():
        if len(values["Cluster"]) == 0:
            values["Cluster"].append("-")


def specificity(contig_dict, table, tag):
    for line in table:
        contig = line.strip()
        contig_dict[contig]["{tag}_specific".format(tag=tag)].append("*")

    for contig, values in contig_dict.items():
        if len(values["{tag}_specific".format(tag=tag)]) == 0:
            values["{tag}_specific".format(tag=tag)].append("-")


def expression(contig_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        contig, red, cer, mar = description[0], float(description[1]), float(description[2]), float(description[3])
        contig_dict[contig]["Rediae_scaledTPM"] += red
        contig_dict[contig]["Cercariae_scaledTPM"] += cer
        contig_dict[contig]["Marita_scaledTPM"] += mar


def write_output_files(contig_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("### Contig_ID - ID of assembled and selected contig\n"
                          "### Orthogroup - ID of the orthogroup group in which the sequence was included\n"
                          "### NCBInt - best hit with NCBI nucleotide database\n"
                          "### NCBInr - best hit with NCBI non-redundant protein database\n"
                          "### SwissProt - best hit with SwissProt protein database\n"
                          "### NCBInt|NCBInr|SwissProt_identity - best hit identity percentage\n"
                          "### Domain_arch - reconstructed protein architecture based on analysis "
                          " of comparison results with a PfamA database\n"
                          "### EggNOG:Preferred_name - Predicted gene name\n"
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
                          "### Rediae|Cercariae|Marita_scaledTPM - the value of expression in rediae, cercariae "
                          "or adult worm stages, respectively, in TPM (Transcripts Per Million) units "
                          "which were between-sample normalized and averaged between two biological replicates\n"
                          "### Rediae|Cercariae|Marita-specific - '*' means that the sequence was attributed to rediae,"
                          "cercariae or adult worm stage specific gene set, respectively, according to Jongeneel`s "
                          "specificity measure\n"
                          "### Cluster - co-expressed gene cluster identifier\n")
        output_file.write("Contig_ID\tOrthogroup\tEggNOG:Preferred_name\tNCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\tDomains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tRediae_scaledTPM\tCercariae_scaledTPM\t"
                          "Marita_scaledTPM\tRediae-specific\tCercariae-specific\t"
                          "Marita-specific\tCluster\n")
        for contig, values in contig_dict.items():
            output_file.write("{contig}\t{ortho}\t{name}\t{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t{swiss}\t"
                              "{swiss_identity}\t{domains}\t{go}\t{ec}\t{ko}\t{kegg}\t{module}\t{reaction}\t"
                              "{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\t"
                              "{red_tpm}\t{cer_tpm}\t{mar_tpm}\t{red_spec}\t{cer_spec}\t{mar_spec}\t"
                              "{cluster}\n".format(contig=contig,
                               ortho=values["Ortho"][0], name=values["EggNOG:Preferred_name"][0],
                               nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                               nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                               swiss=values["Anno_swiss"]["hit"][0], swiss_identity=values["Anno_swiss"]["identity"][0],
                               domains="|".join(values["domains"]),
                               go=values["EggNOG:GO_terms"][0], ec=values["EggNOG:EC_number"][0],
                               ko=values["EggNOG:KEGG_KO"][0], kegg=values["EggNOG:KEGG_Pathway"][0],
                               module=values["EggNOG:KEGG_Module"][0], reaction=values["EggNOG:KEGG_Reaction"][0],
                               rclass=values["EggNOG:rclass"][0], brite=values["EggNOG:BRITE"][0],
                               tc=values["EggNOG:KEGG_TC"][0], cazy=values["EggNOG:CAZy"][0],
                               bigg=values["EggNOG:BiGG_Reaction"][0], og=values["EggNOG:OG"][0],
                               cog=values["EggNOG:COG_cat"][0], description=values["EggNOG:Description"][0],
                               red_tpm=values["Rediae_scaledTPM"], cer_tpm=values["Cercariae_scaledTPM"],
                               mar_tpm=values["Marita_scaledTPM"], red_spec=values["Rediae_specific"][0],
                               cer_spec=values["Cercaria_specific"][0], mar_spec=values["Marita_specific"][0],
                               cluster=values["Cluster"][0]))


if __name__ == "__main__":
    contig_dict, clusters_dict = {}, {}
    print("***** Fasta files parsing *****")
    nucl_parsing(contig_dict, args.nucl)
    amino_parsing(contig_dict, args.amino)
    print("***** Orthogroups parsing *****")
    orthogroups_parsing(contig_dict, args.ortho)
    print("***** eggNOG-mapper output parsing *****")
    eggNOG_mapper_parsing(contig_dict, args.eggnog)
    print("***** BLAST results parsing *****")
    BLAST_annotation(contig_dict, args.blast_swiss, "swiss")
    BLAST_annotation(contig_dict, args.blast_nt, "nt")
    BLAST_annotation(contig_dict, args.blast_nr, "nr")
    print("***** Domain architecture parsing *****")
    domains_parsing(contig_dict, args.domains)
    print("***** Cluster parsing *****")
    clusters(contig_dict, args.clusters, clusters_dict)
    print("***** Table with averaged expression values parsing *****")
    expression(contig_dict, args.scaled_tpm)
    print("***** Files with IDs of stage-specific sequences parsing *****")
    specificity(contig_dict, args.rediae_specific, "Rediae")
    specificity(contig_dict, args.cercariae_specific, "Cercaria")
    specificity(contig_dict, args.marita_specific, "Marita")
    print("***** Output file creating *****")
    write_output_files(contig_dict, args.output)
