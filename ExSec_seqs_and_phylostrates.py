try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--classic_exsec', type=argparse.FileType('r'), required=True,
                    help="Fasta file with potential classical excretory/secretory sequences")
parser.add_argument('--nonclassic_exsec', type=argparse.FileType('r'), required=True,
                    help="Fasta file with potential nonclassical excretory/secretory sequences")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[protein_ID] = mrca_name


def fasta_parsing(fasta, fasta_list):
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for fasta in fasta_seqs:
        fasta_list.append(fasta.id)


def summarization(phylostratr_dict, classic_list, nonclassic_list, summary_dict):
    phylostrates = set([value for value in phylostratr_dict.values()])
    for phylostrate in phylostrates:
        summary_dict[phylostrate] = {"classic": [], "nonclassic": []}

    for classic_seq in classic_list:
        summary_dict[phylostratr_dict[classic_seq]]["classic"].append(classic_seq)

    for nonclassic_seq in nonclassic_list:
        summary_dict[phylostratr_dict[nonclassic_seq]]["nonclassic"].append(nonclassic_seq)


def output_writing(output, classic_list, nonclassic_list, summary_dict):
    with open("{output}.exsec_seqs_in_phylostrates_summary.protein_IDs.tsv".format(
            output=output), 'a') as summary_with_IDs:
        summary_with_IDs.write("Phylostrates\tClassical_ExSec\tNonclassical_ExSec\n")
        for phylostrate, values in summary_dict.items():
            summary_with_IDs.write("{phylo}\t{classic}\t{nonclassic}\n".format(
                phylo=phylostrate, classic=";".join(values["classic"]), nonclassic=";".join(values["nonclassic"])))

    with open("{output}.exsec_seqs_in_phylostrates_summary.protein_counts_and_percents.tsv".format(
            output=output), 'a') as summary_with_counts:
        summary_with_counts.write("Phylostrates\tClassical_ExSec\tNonclassical_ExSec\n")
        for phylostrate, values in summary_dict.items():
            summary_with_counts.write("{phylo}\t{classic}\t{nonclassic}\n".format(
                phylo=phylostrate,
                classic="{count} ({percent}%)".format(
                    count=len(values["classic"]),
                    percent=round((len(values["classic"])/len(classic_list))*100, 2)),
                nonclassic="{count} ({percent}%)".format(
                    count=len(values["nonclassic"]),
                    percent=round((len(values["nonclassic"])/len(nonclassic_list))*100, 2))
            ))


if __name__ == "__main__":
    phylostratr_dict, summary_dict, classic_list, nonclassic_list = {}, {}, [], []
    print("***** Input files parsing *****")
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    fasta_parsing(args.classic_exsec, classic_list)
    fasta_parsing(args.nonclassic_exsec, nonclassic_list)
    print("***** Data analysis and summarization *****")
    summarization(phylostratr_dict, classic_list, nonclassic_list, summary_dict)
    print("***** Output file writing *****")
    output_writing(args.output, classic_list, nonclassic_list, summary_dict)