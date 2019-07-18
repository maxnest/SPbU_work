try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

try:
    import numpy
except ImportError:
    print("Please check if module 'numpy' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast', type=argparse.FileType('r'), required=True)
parser.add_argument('--KOBAS_Hsapiens', type=argparse.FileType('r'), required=True)
parser.add_argument('--GO', type=argparse.FileType('r'), required=True)
parser.add_argument('--sites', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for different sites")
parser.add_argument('--head_4', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for HEAD in 4 hours")
parser.add_argument('--head_12', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for HEAD in 12 hours")
parser.add_argument('--head_24', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for HEAD in 24 hours")
parser.add_argument('--head_48', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for HEAD in 48 hours")
parser.add_argument('--head_96', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for HEAD in 96 hours")
parser.add_argument('--tail_4', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for TAIL in 4 hours")
parser.add_argument('--tail_12', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for TAIL in 12 hours")
parser.add_argument('--tail_24', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for TAIL in 24 hours")
parser.add_argument('--tail_48', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for TAIL in 48 hours")
parser.add_argument('--tail_96', type=argparse.FileType('r'), required=True,
                    help="DESeq2 table with LFC for TAIL in 96 hours")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ID": [seq.id.split(" ")[0]], "nucl": seq.seq, "amino": "",
                               "domains": [], "GO_bio": [], "GO_mol": [], "GO_cell": [], "KEGG_Hsapiens": [],
                               "Anno_swiss": {"hit": [], "identity": []},
                               "Tail_vs_Head": [], "Head_0_vs_4h": [], "Head_0_vs_12h": [], "Head_0_vs_24h": [],
                               "Head_0_vs_48h": [], "Head_0_vs_96h": [], "Tail_0_vs_4h": [], "Tail_0_vs_12h": [],
                               "Tail_0_vs_24h": [], "Tail_0_vs_48h": [], "Tail_0_vs_96h": [], "baseMean": [],
                               "padj": []}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".p")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"


def GO_parsing(contig_dict, GO):
    head = GO.readline()
    for line in GO:
        description = line.strip().split("\t")
        contig_ID, category, ID, rank, anno = description[0].strip().split(".p")[0], description[1], \
                                              "GO:{ID}".format(ID=description[2]), description[-1], description[3]
        if contig_ID in contig_dict.keys():
            if rank == "1":
                if category == "MF":
                    contig_dict[contig_ID]["GO_mol"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
                elif category == "CC":
                    contig_dict[contig_ID]["GO_cell"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
                elif category == "BP":
                    contig_dict[contig_ID]["GO_bio"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))

    for contig, values in contig_dict.items():
        if len(values["GO_mol"]) == 0:
            contig_dict[contig]["GO_mol"].append("-")

        if len(values["GO_cell"]) == 0:
            contig_dict[contig]["GO_cell"].append("-")

        if len(values["GO_bio"]) == 0:
            contig_dict[contig]["GO_bio"].append("-")


def BLAST_annotation(contig_dict, BLAST, key_tag):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        ID, hit_name = description[0].strip().split(".p")[0], description[3]
        if ID in contig_dict.keys():
            if hit_name == "no hits":
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
            else:
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append(description[8])

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


def KEGG_parsing(contig_dict, KEGG, key_tag):
    for number in range(5):
        head = KEGG.readline()

    for line in KEGG:
        description = line.strip().split("\t")
        if len(description) > 1:
            ID, hit = description[0].strip().split(".p")[0], description[1]
            if ID in contig_dict.keys():
                if hit == "None":
                    contig_dict[ID][key_tag].append("No hit")
                else:
                    KEGG_hit = hit.split("||")
                    contig_dict[ID][key_tag].append(KEGG_hit[0])
        else:
            break

    for contig in contig_dict.keys():
        if len(contig_dict[contig][key_tag]) == 0:
            contig_dict[contig][key_tag].append("-")


#   def deseq2_matrix(contig_dict, matrix):
#       head = matrix.readline()
#       for line in matrix:
#            description = line.strip().split(",")
#            ID, intercept, sites, HEAD_4, HEAD_12, HEAD_24, HEAD_48, HEAD_96, TAIL_4, TAIL_12, TAIL_24,
#            TAIL_48, TAIL_96 = description[0].split('"')[1], description[1], description[2], description[3],
#            description[4], description[5], description[6], description[7], description[8], description[9],
#            description[10], description[11], description[12]
#            contig_dict[ID]["Intercept"].append(intercept)
#            contig_dict[ID]["Tail_vs_Head"].append(sites)
#            contig_dict[ID]["Head_0_vs_4h"].append(HEAD_4)
#            contig_dict[ID]["Head_0_vs_12h"].append(HEAD_12)
#            contig_dict[ID]["Head_0_vs_24h"].append(HEAD_24)
#            contig_dict[ID]["Head_0_vs_48h"].append(HEAD_48)
#            contig_dict[ID]["Head_0_vs_96h"].append(HEAD_96)
#            contig_dict[ID]["Tail_0_vs_4h"].append(TAIL_4)
#            contig_dict[ID]["Tail_0_vs_12h"].append(TAIL_12)
#            contig_dict[ID]["Tail_0_vs_24h"].append(TAIL_24)
#            contig_dict[ID]["Tail_0_vs_48h"].append(TAIL_48)
#            contig_dict[ID]["Tail_0_vs_96h"].append(TAIL_96)

#       for contig in contig_dict.keys():
#           for el in ["Intercept", "Tail_vs_Head", "Head_0_vs_4h", "Head_0_vs_12h", "Head_0_vs_24h", "Head_0_vs_48h",
#                      "Head_0_vs_96h", "Tail_0_vs_4h", "Tail_0_vs_12h", "Tail_0_vs_24h", "Tail_0_vs_48h",
#                      "Tail_0_vs_96h"]:
#               if len(contig_dict[contig][el]) == 0:
#                   contig_dict[contig][el].append("low")


def deseq2_results(contig_dict, results, tag):
    head = results.readline()
    for line in results:
        description = line.strip().split(",")
        ID, baseMean, LFC, padj = description[0].split('"')[1], description[1], description[2], description[-1]
        if len(contig_dict[ID]["baseMean"]) == 0:
            contig_dict[ID]["baseMean"].append(baseMean)

        if len(contig_dict[ID]["padj"]) == 0:
            contig_dict[ID]["padj"].append(padj)

        if len(contig_dict[ID][tag]) == 0:
            contig_dict[ID][tag].append(LFC)

    for contig in contig_dict.keys():
        for el in ["baseMean", tag, "padj"]:
            if len(contig_dict[contig][el]) == 0:
                contig_dict[contig][el].append("low")


def write_output(contig_dict, output_tag):
    with open("{output_tag}.tab".format(output_tag=output_tag), 'a') as output:
        output.write("Contig_ID\tSwiss\tIdentity\tKEGG_Hsapiens\tGO_bio\tGO_mol\tGO_cell\tbaseMean\tpadj\t"
                     "Tail_vs_Head\tHead_0_vs_4h\tHead_0_vs_12h\tHead_0_vs_24h\tHead_0_vs_48h\tHead_0_vs_96h\t"
                     "Tail_0_vs_4h\tTail_0_vs_12h\tTail_0_vs_24h\tTail_0_vs_48h\tTail_0_vs_96h\tNucl_seq\tAA_seq"
                     "\tDomains\n")

        for contig, values in contig_dict.items():
            output.write("{ID}\t{Swiss}\t{Identity}\t{KEGG}\t{bio}\t{mol}\t{cell}\t{base}\t{padj}\t"
                         "{Sites}\t{Head_4h}\t{Head_12h}\t{Head_24h}\t{Head_48h}\t{Head_96h}\t{Tail_4h}\t{Tail_12h}\t"
                         "{Tail_24h}\t{Tail_48h}\t{Tail_96h}\t{Nucl}\t{AA}\t{Domains}\n".format(ID=contig,
            Swiss=values["Anno_swiss"]["hit"][0], Identity=values["Anno_swiss"]["identity"][0],
            KEGG=values["KEGG_Hsapiens"][0], bio=values["GO_bio"][0], mol=values["GO_mol"][0],
            cell=values["GO_cell"][0], base=values["baseMean"][0], padj=values["padj"][0],
            Sites=values["Tail_vs_Head"][0], Head_4h=values["Head_0_vs_4h"][0],  Head_12h=values["Head_0_vs_12h"][0],
            Head_24h=values["Head_0_vs_24h"][0], Head_48h=values["Head_0_vs_48h"][0], Head_96h=values["Head_0_vs_96h"][0],
            Tail_4h=values["Tail_0_vs_4h"][0], Tail_12h=values["Tail_0_vs_12h"][0], Tail_24h=values["Tail_0_vs_24h"][0],
            Tail_48h=values["Tail_0_vs_48h"][0], Tail_96h=values["Tail_0_vs_96h"][0], Nucl=values["nucl"],
            AA=values["amino"], Domains="|".join(values["domains"])))


if __name__ == "__main__":
    contig_dict = {}
    print("Parsing fasta-file with nucleotide sequences")
    nucl_parsing(contig_dict, args.nucl)
    print("Parsing fasta-file with aminoacid sequences")
    amino_parsing(contig_dict, args.amino)
    print("Parsing table with Gene Ontology")
    GO_parsing(contig_dict, args.GO)
    print("Parsing table with blast results")
    BLAST_annotation(contig_dict, args.blast, "swiss")
    print("Parsing table with domains")
    domains_parsing(contig_dict, args.domains)
    print("Parsing table with pathways")
    KEGG_parsing(contig_dict, args.KOBAS_Hsapiens, "KEGG_Hsapiens")
    print("Parsing tables with DE results")
    deseq2_results(contig_dict, args.sites, "Tail_vs_Head")
    deseq2_results(contig_dict, args.head_4, "Head_0_vs_4h")
    deseq2_results(contig_dict, args.head_12, "Head_0_vs_12h")
    deseq2_results(contig_dict, args.head_24, "Head_0_vs_24h")
    deseq2_results(contig_dict, args.head_48, "Head_0_vs_48h")
    deseq2_results(contig_dict, args.head_96, "Head_0_vs_96h")
    deseq2_results(contig_dict, args.tail_4, "Tail_0_vs_4h")
    deseq2_results(contig_dict, args.tail_12, "Tail_0_vs_12h")
    deseq2_results(contig_dict, args.tail_24, "Tail_0_vs_24h")
    deseq2_results(contig_dict, args.tail_48, "Tail_0_vs_48h")
    deseq2_results(contig_dict, args.tail_96, "Tail_0_vs_96h")
    print("Output file creating")
    write_output(contig_dict, args.output)
