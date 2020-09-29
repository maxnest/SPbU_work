try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--classical_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--nonclassical_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum phylostratigraphy analysis")
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="Table with description of phylostratigraphic levels")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def levels_parsing(levels, levels_dict):
    for line in levels:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def phylostratr_table_parsing(phylostratr, levels_dict, phylostrata_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostrata_dict[protein_ID] = levels_dict[mrca_name]


def fasta_parser(fasta, fasta_list):
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for fasta in fasta_seqs:
        fasta_list.append(fasta.id)


def results_merging(phylostrata_dict, classical_exsec_list, nonclassical_exsec_list, merged_dict):
    for protein, level in phylostrata_dict.items():
        if level not in merged_dict.keys():
            merged_dict[level] = {"classical": [], "nonclassical": []}

    for protein in classical_exsec_list:
        merged_dict[phylostrata_dict[protein]]["classical"].append(protein)

    for protein in nonclassical_exsec_list:
        merged_dict[phylostrata_dict[protein]]["nonclassical"].append(protein)


def output_writing(output, classical_exsec_list, nonclassical_exsec_list, merged_dict):
    with open("{output}_exsec_seq_in_phylostrates.tsv".format(output=output), 'a') as output_file:
        output_file.write("Phylostratas\tClassical_exsec:count\tClassical_exsec:percent\t"
                          "Nonclassical_exsec:count\tNonclassical_exsec:percent\n")
        for level, values in merged_dict.items():
            output_file.write("{level}\t{classical_count}\t{classical_percent}\t"
                              "{nonclassical_count}\t{nonclassical_percent}\n".format(
                                level=level,
                                classical_count=len(values["classical"]),
                                classical_percent=round((len(values["classical"])/len(classical_exsec_list))*100, 2),
                                nonclassical_count=len(values["nonclassical"]),
                                nonclassical_percent=
                                round((len(values["nonclassical"])/len(nonclassical_exsec_list))*100, 2)))


if __name__ == "__main__":
    levels_dict, phylostrata_dict, merged_dict = {}, {}, {}
    classical_exsec_list, nonclassical_exsec_list = [], []
    print("***** Input files parsing *****")
    levels_parsing(args.levels, levels_dict)
    phylostratr_table_parsing(args.phylostratr, levels_dict, phylostrata_dict)
    fasta_parser(args.classical_exsec, classical_exsec_list)
    fasta_parser(args.nonclassical_exsec, nonclassical_exsec_list)
    print("***** Results merging *****")
    results_merging(phylostrata_dict, classical_exsec_list, nonclassical_exsec_list, merged_dict)
    print("***** Output writing *****")
    output_writing(args.output, classical_exsec_list, nonclassical_exsec_list, merged_dict)
