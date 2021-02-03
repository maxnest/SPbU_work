try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np
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
parser.add_argument('--animaltf', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--averaged_expression', type=argparse.FileType('r'), required=True)
parser.add_argument('--RNentropy_over', type=argparse.FileType('r'), required=True)
parser.add_argument('--RNentropy_under', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_ts', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_ts', type=argparse.FileType('r'), required=True)
parser.add_argument('--both_sites_ts', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True)
parser.add_argument('--rbbh', type=argparse.FileType('r'), required=True)
parser.add_argument('--second_species_tag', type=str, required=True)
parser.add_argument('--nucl_fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino_fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def gene_map_parsing(gene_dict, gene_map):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_ID, transcript_ID, protein_ID = description[0], description[1], description[2]
        gene_dict[gene_ID] = {"transcript": transcript_ID, "protein": protein_ID,
                              "domains": [],
                              "Anno_swiss": {"hit": [], "identity": []},
                              "Anno_nt": {"hit": [], "identity": []},
                              "Anno_nr": {"hit": [], "identity": []},
                              "Anno_neuropep": {"hit": [], "identity": []},
                              "Anno_animalTF": {"hit": [], "identity": []},
                              "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                              "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                              "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                              "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                              "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                              "Head_0h_TPMs": 0, "Head_4h_TPMs": 0, "Head_12h_TPMs": 0,
                              "Head_24h_TPMs": 0, "Head_48h_TPMs": 0, "Head_96h_TPMs": 0,
                              "Tail_0h_TPMs": 0, "Tail_4h_TPMs": 0, "Tail_12h_TPMs": 0,
                              "Tail_24h_TPMs": 0, "Tail_48h_TPMs": 0, "Tail_96h_TPMs": 0,
                              "Head_cluster": [], "Tail_cluster": [],
                              "Head_TS": [], "Tail_TS": [], "Both_sites_TS": [],
                              "Phylostrata": [], "RBBH": [],
                              "RNentropy_over": [], "RNentropy_under": [], "Ortho": [],
                              "Nucl_seq": "", "Amino_seq": ""}


def fasta_parsing(fasta_file, tag, seq_2_gene_dict, gene_dict):
    fasta = SeqIO.parse(fasta_file, 'fasta')
    for seq in fasta:
        gene_dict[seq_2_gene_dict[seq.id]]["{tag}_seq".format(tag=tag)] += seq.seq


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
                values[key].append("-")


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


def RNentropy_parsing(gene_dict, table, tag):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, sample = description[0][1:-1], description[-1][1:-1]
        gene_dict[gene]["RNentropy_{tag}".format(tag=tag)].append(sample)

    for gene, values in gene_dict.items():
        if len(values["RNentropy_{tag}".format(tag=tag)]) == 0:
            values["RNentropy_{tag}".format(tag=tag)].append("-")


def expression_parsing(gene_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, head_0h, head_4h,	 head_12h, head_24h, head_48h, head_96h, \
        tail_0h, tail_4h, tail_12h, tail_24h, tail_48h, tail_96h = description[0], \
            float(description[1]), float(description[2]), float(description[3]), float(description[4]), \
            float(description[5]), float(description[6]), float(description[7]), float(description[8]), \
            float(description[9]), float(description[10]), float(description[11]), float(description[12])
        gene_dict[gene]["Head_0h_TPMs"] += head_0h
        gene_dict[gene]["Head_4h_TPMs"] += head_4h
        gene_dict[gene]["Head_12h_TPMs"] += head_12h
        gene_dict[gene]["Head_24h_TPMs"] += head_24h
        gene_dict[gene]["Head_48h_TPMs"] += head_48h
        gene_dict[gene]["Head_96h_TPMs"] += head_96h
        gene_dict[gene]["Tail_0h_TPMs"] += tail_0h
        gene_dict[gene]["Tail_4h_TPMs"] += tail_4h
        gene_dict[gene]["Tail_12h_TPMs"] += tail_12h
        gene_dict[gene]["Tail_24h_TPMs"] += tail_24h
        gene_dict[gene]["Tail_48h_TPMs"] += tail_48h
        gene_dict[gene]["Tail_96h_TPMs"] += tail_96h


def clusters(gene_dict, table, tag):
    header = table.readline().strip().split("\t")
    cluster_dict = {}

    for el in header:
        cluster_dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                cluster_dict[header[genes.index(gene)]].append(gene)

    for gene in gene_dict.keys():
        for cluster, genes in cluster_dict.items():
            if gene in genes:
                gene_dict[gene]["{tag}_cluster".format(tag=tag)].append(cluster)

    for gene, values in gene_dict.items():
        if len(values["{tag}_cluster".format(tag=tag)]) == 0:
            values["{tag}_cluster".format(tag=tag)].append("-")


def time_series_parsing(gene_dict, ts, tag):
    for line in ts:
        gene = line.strip()
        gene_dict[gene][tag].append("*")

    for gene, values in gene_dict.items():
        if len(values[tag]) == 0:
            values[tag].append("-")


def phylostratr_parsing(gene_dict, prot_2_gene_dict, phylostratr):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        if protein_ID in prot_2_gene_dict.keys():
            gene_dict[prot_2_gene_dict[protein_ID]]["Phylostrata"].append("{ps}:{mrca_name}".format(
                ps=ps, mrca_name=mrca_name))


def RBBH_parsing(gene_dict, rbbh):
    header = rbbh.readline()
    for line in rbbh:
        description = line.strip().split("\t")
        rbbh_ID, first_gene, second_gene = description[0], description[1], description[2]
        if first_gene in gene_dict.keys():
            gene_dict[first_gene]["RBBH"].append("{id}:{second}".format(id=rbbh_ID, second=second_gene))

        if second_gene in gene_dict.keys():
            gene_dict[second_gene]["RBBH"].append("{id}:{first}".format(id=rbbh_ID, first=first_gene))

    for gene, values in gene_dict.items():
        if len(values["RBBH"]) == 0:
            values["RBBH"].append("-")


def write_output_files(gene_dict, output, second_species_tag):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tTranscript_ID\tProtein_ID\tPhylostrates\tOrthogroup\tEggNOG:Preferred_name\t"
                          "NCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\tNeuroPep\tNeuroPep_identity\t"
                          "AnimalTF\tAnimalTF_identity\tDomains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tRBBH_to_{tag}\t"
                          "Head_0h_log2_averaged_TPMs\tHead_4h_log2_averaged_TPMs\t"
                          "Head_12h_log2_averaged_TPMs\tHead_24h_log2_averaged_TPMs\t"
                          "Head_48h_log2_averaged_TPMs\tHead_96h_log2_averaged_TPMs\t"
                          "Tail_0h_log2_averaged_TPMs\tTail_4h_log2_averaged_TPMs\t"
                          "Tail_12h_log2_averaged_TPMs\tTail_24h_log2_averaged_TPMs\t"
                          "Tail_48h_log2_averaged_TPMs\tTail_96h_log2_averaged_TPMs\t"
                          "Overexpression_in_samples\tUnderexpression_in_samples\t"
                          "DiffExp_in_time_series:head\tDiffExp_in_time_series:tail\t"
                          "DiffExp_in_time_series:between_sites\t"
                          "Co-expression_clusters:head\tCo-expression_clusters:tail\t"
                          "Nucl_seqs\tAmino_seqs\n".format(tag=second_species_tag))
        for gene, values in gene_dict.items():
            output_file.write("{gene}\t{trans}\t{prot}\t{phylostrates}\t{ortho}\t{name}\t"
                              "{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t"
                              "{swiss}\t{swiss_identity}\t{neuro}\t{neuro_identity}\t"
                              "{tf}\t{tf_identity}\t{domains}\t{go}\t{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t"
                              "{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\t{rbbh}\t"
                              "{head_0h}\t{head_4h}\t{head_12h}\t{head_24h}\t{head_48h}\t{head_96h}\t"
                              "{tail_0h}\t{tail_4h}\t{tail_12h}\t{tail_24h}\t{tail_48h}\t{tail_96h}\t"
                              "{overexp}\t{underexp}\t{ts_head}\t{ts_tail}\t{ts_sites}\t"
                              "{head_clusters}\t{tail_clusters}\t{nucl}\t{amino}\n".format(
                                gene=gene, trans=values["transcript"], prot=values["protein"],
                                phylostrates=values["Phylostrata"][0], ortho=values["Ortho"][0],
                                name=values["EggNOG:Preferred_name"][0],
                                nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                                nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                                swiss=values["Anno_swiss"]["hit"][0], swiss_identity=values["Anno_swiss"]["identity"][0],
                                neuro=values["Anno_neuropep"]["hit"][0],
                                neuro_identity=values["Anno_neuropep"]["identity"][0],
                                tf=values["Anno_animalTF"]["hit"][0], tf_identity=values["Anno_animalTF"]["identity"][0],
                                domains="|".join(values["domains"]), go=values["EggNOG:GO_terms"][0],
                                ec=values["EggNOG:EC_number"][0], ko=values["EggNOG:KEGG_KO"][0],
                                pathway=values["EggNOG:KEGG_Pathway"][0], module=values["EggNOG:KEGG_Module"][0],
                                reaction=values["EggNOG:KEGG_Reaction"][0], rclass=values["EggNOG:rclass"][0],
                                brite=values["EggNOG:BRITE"][0], tc=values["EggNOG:KEGG_TC"][0],
                                cazy=values["EggNOG:CAZy"][0], bigg=values["EggNOG:BiGG_Reaction"][0],
                                og=values["EggNOG:OG"][0], cog=values["EggNOG:COG_cat"][0],
                                description=values["EggNOG:Description"][0], rbbh=values["RBBH"][0],
                                head_0h=np.log2(values["Head_0h_TPMs"] + 1),
                                head_4h=np.log2(values["Head_4h_TPMs"] + 1),
                                head_12h=np.log2(values["Head_12h_TPMs"] + 1),
                                head_24h=np.log2(values["Head_24h_TPMs"] + 1),
                                head_48h=np.log2(values["Head_48h_TPMs"] + 1),
                                head_96h=np.log2(values["Head_96h_TPMs"] + 1),
                                tail_0h=np.log2(values["Tail_0h_TPMs"] + 1),
                                tail_4h=np.log2(values["Tail_4h_TPMs"] + 1),
                                tail_12h=np.log2(values["Tail_12h_TPMs"] + 1),
                                tail_24h=np.log2(values["Tail_24h_TPMs"] + 1),
                                tail_48h=np.log2(values["Tail_48h_TPMs"] + 1),
                                tail_96h=np.log2(values["Tail_96h_TPMs"] + 1),
                                overexp=";".join(values["RNentropy_over"]),
                                underexp=";".join(values["RNentropy_under"]),
                                ts_head=values["Head_TS"][0],
                                ts_tail=values["Tail_TS"][0], ts_sites=values["Both_sites_TS"][0],
                                head_clusters=values["Head_cluster"][0], tail_clusters=values["Tail_cluster"][0],
                                nucl=values["Nucl_seq"], amino=values["Amino_seq"]))


if __name__ == "__main__":
    gene_dict = {}
    gene_map_parsing(gene_dict, args.gene_map)
    prot_2_gene_dict = {gene_dict[gene]["protein"]: gene for gene in gene_dict.keys()}
    trans_2_gene_dict = {gene_dict[gene]["transcript"]: gene for gene in gene_dict.keys()}
    print("***** Fasta files parsing *****")
    fasta_parsing(args.nucl_fasta, "Nucl", trans_2_gene_dict, gene_dict)
    fasta_parsing(args.amino_fasta, "Amino", prot_2_gene_dict, gene_dict)
    print("***** Orthogroups parsing *****")
    orthogroups_parsing(gene_dict, prot_2_gene_dict, args.ortho)
    print("***** eggNOG-mapper output parsing *****")
    eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, args.eggnog)
    print("***** BLAST results parsing *****")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_swiss, "swiss")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nt, "nt")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nr, "nr")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.neuropep, "neuropep")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.animaltf, "animalTF")
    print("***** Domain architecture parsing *****")
    domains_parsing(gene_dict, prot_2_gene_dict, args.domains)
    print("***** Expression parsing *****")
    expression_parsing(gene_dict, args.averaged_expression)
    RNentropy_parsing(gene_dict, args.RNentropy_over, "over")
    RNentropy_parsing(gene_dict, args.RNentropy_under, "under")
    time_series_parsing(gene_dict, args.head_ts, "Head_TS")
    time_series_parsing(gene_dict, args.tail_ts, "Tail_TS")
    time_series_parsing(gene_dict, args.both_sites_ts, "Both_sites_TS")
    clusters(gene_dict, args.head_clusters, "Head")
    clusters(gene_dict, args.tail_clusters, "Tail")
    print("***** Phylostratigraphic analysis results parsing *****")
    phylostratr_parsing(gene_dict, prot_2_gene_dict, args.phylostratr)
    print("***** RBBHs parsing *****")
    RBBH_parsing(gene_dict, args.rbbh)
    print("***** Output writing *****")
    write_output_files(gene_dict, args.output, args.second_species_tag)