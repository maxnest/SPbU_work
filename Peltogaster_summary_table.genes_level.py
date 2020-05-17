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
parser.add_argument('--neuropep', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--scaled_tpm', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_body_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_markers', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_markers', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_markers', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_markers', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_body_markers', type=argparse.FileType('r'), required=True)
parser.add_argument('--classic_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--nonclassic_exsec', type=argparse.FileType('r'), required=True)
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
                              "Anno_nr": {"hit": [], "identity": []}, "Anno_neuropep": {"hit": [], "identity": []},
                              "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                              "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                              "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                              "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                              "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                              "Externa_scaledTPM": 0, "Growing_scaledTPM": 0, "Middle_scaledTPM": 0,
                              "Terminal_scaledTPM": 0, "Whole_body_scaledTPM": 0,
                              "Specificity": [], "Markers": [], "Ortho": [], "ExSec": []}


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
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["hit"].append("-")
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["identity"].append("-")
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


def specificity(gene_dict, table, tag):
    for line in table:
        gene = line.strip()
        gene_dict[gene]["Specificity"].append(tag)


def markers(gene_dict, table, tag):
    for line in table:
        gene = line.strip()
        gene_dict[gene]["Markers"].append(tag)


def expression(gene_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, externa, growing, middle, terminal, whole_body = \
            description[0], float(description[1]), float(description[2]), float(description[3]), \
            float(description[4]), float(description[5])
        gene_dict[gene]["Externa_scaledTPM"] += externa
        gene_dict[gene]["Growing_scaledTPM"] += growing
        gene_dict[gene]["Middle_scaledTPM"] += middle
        gene_dict[gene]["Terminal_scaledTPM"] += terminal
        gene_dict[gene]["Whole_body_scaledTPM"] += whole_body


def exsec_parsing(gene_dict, exsec_seqs, prot_2_gene_dict, tag):
    for fasta in exsec_seqs:
        description = fasta.id.strip().split(" ")
        gene_dict[prot_2_gene_dict[description[0]]]["ExSec"].append(tag)


def write_output_files(gene_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tTranscript_ID\tProtein_ID\tOrthogroup\tEggNOG:Preferred_name\t"
                          "NCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\tNeuroPep\tNeuroPep_identity\t"
                          "Domains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tExterna_scaledTPM\t"
                          "Growing_trunk_scaledTPM\tMain_trunk_part_scaledTPM\tThoracic_part_of_interna\t"
                          "Whole_body_scaledTPM\tSpecificity\tMolecular_markers\tSecretory\Excretory\n")
        for gene, values in gene_dict.items():
            output_file.write("{gene}\t{trans}\t{prot}\t{ortho}\t{name}\t{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t"
                              "{swiss}\t{swiss_identity}\t{neuro}\t{neuro_identity}\t{domains}\t{go}\t"
                              "{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t"
                              "{og}\t{cog}\t{description}\t{externa}\t{growing}\t{middle}\t{terminal}\t{whole}\t"
                              "{specificity}\t{markers}\t{exsec}\n".format(
                                gene=gene, trans=values["transcript"], prot=values["protein"],
                                ortho=values["Ortho"][0], name=values["EggNOG:Preferred_name"][0],
                                nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                                nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                                swiss=values["Anno_swiss"]["hit"][0],
                                swiss_identity=values["Anno_swiss"]["identity"][0],
                                neuro=values["Anno_neuropep"]["hit"][0],
                                neuro_identity=values["Anno_neuropep"]["identity"][0],
                                domains="|".join(values["domains"]), go=values["EggNOG:GO_terms"][0],
                                ec=values["EggNOG:EC_number"][0], ko=values["EggNOG:KEGG_KO"][0],
                                pathway=values["EggNOG:KEGG_Pathway"][0], module=values["EggNOG:KEGG_Module"][0],
                                reaction=values["EggNOG:KEGG_Reaction"][0], rclass=values["EggNOG:rclass"][0],
                                brite=values["EggNOG:BRITE"][0], tc=values["EggNOG:KEGG_TC"][0],
                                cazy=values["EggNOG:CAZy"][0], bigg=values["EggNOG:BiGG_Reaction"][0],
                                og=values["EggNOG:OG"][0], cog=values["EggNOG:COG_cat"][0],
                                description=values["EggNOG:Description"][0], externa=values["Externa_scaledTPM"],
                                growing=values["Growing_scaledTPM"], middle=values["Middle_scaledTPM"],
                                terminal=values["Terminal_scaledTPM"], whole=values["Whole_body_scaledTPM"],
                                specificity=values["Specificity"][0], markers=values["Markers"][0],
                                exsec=values["ExSec"][0]))


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
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.neuropep, "neuropep")
    print("***** Domain architecture parsing *****")
    domains_parsing(gene_dict, prot_2_gene_dict, args.domains)
    print("***** Files with IDs of sequences with specific expression patterns parsing *****")
    specificity(gene_dict, args.externa_specific, "Externa")
    specificity(gene_dict, args.growing_specific, "Growing_trunk")
    specificity(gene_dict, args.middle_specific, "Main_trunk_part")
    specificity(gene_dict, args.terminal_specific, "Thoracic_part_of_interna")
    specificity(gene_dict, args.whole_body_specific, "Whole_body")
    for gene, values in gene_dict.items():
        if len(values["Specificity"]) == 0:
            values["Specificity"].append('-')
    markers(gene_dict, args.externa_markers, "Externa")
    markers(gene_dict, args.growing_markers, "Growing_trunk")
    markers(gene_dict, args.middle_markers, "Main_trunk_part")
    markers(gene_dict, args.terminal_markers, "Thoracic_part_of_interna")
    markers(gene_dict, args.whole_body_markers, "Whole_body")
    for gene, values in gene_dict.items():
        if len(values["Markers"]) == 0:
            values["Markers"].append('-')
    print("***** Table with averaged expression values parsing *****")
    expression(gene_dict, args.scaled_tpm)
    print("***** Parsing of fasta-files with exsec sequences *****")
    classic_seqs = SeqIO.parse(args.classic_exsec, "fasta")
    nonclassic_seqs = SeqIO.parse(args.nonclassic_exsec, "fasta")
    exsec_parsing(gene_dict, classic_seqs, prot_2_gene_dict, "Classic")
    exsec_parsing(gene_dict, nonclassic_seqs, prot_2_gene_dict, "Nonclassic")
    for gene, values in gene_dict.items():
        if len(values["ExSec"]) == 0:
            values["ExSec"].append('-')
    print("***** Output file creating *****")
    write_output_files(gene_dict, args.output)


