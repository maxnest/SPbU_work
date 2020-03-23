try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--busco_tsv', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def blast_parsing(blast, blast_dict):
    for line in blast:
        description = line.strip().split("\t")
        qseq, sseq, identity = description[0], description[1], float(description[3])
        if qseq != sseq:
            blast_dict["{qseq}_vs_{sseq}".format(qseq=qseq, sseq=sseq)] = identity


def busco_tsv_parsing(busco, busco_dict, busco_pairs, blast_dict, warning_pairs):
    for line in busco:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            if description[1] == "Duplicated":
                ID, sequence = description[0], description[2]
                if ID not in busco_dict.keys():
                    busco_dict[ID] = []
                    busco_dict[ID].append(sequence)
                else:
                    busco_dict[ID].append(sequence)

    for ID, sequences in busco_dict.items():
        pairs = ["_vs_".join(pair) for pair in itertools.combinations(sequences, r=2)]
        busco_pairs[ID] = {pair: blast_dict[pair] for pair in pairs if pair.split("_vs_")[0] != pair.split("_vs_")[1]
                           and pair in blast_dict.keys()}
        for pair in pairs:
            if pair not in blast_dict.keys():
                warning_pairs.append(pair)
    print("***** NB! Script found {num} of warning pairs: *****\n"
          "{warning_pairs}\n"
          "*****".format(num=len(warning_pairs), warning_pairs="\n".join(warning_pairs)))


def summary(busco_pairs, warning_pairs, output):
    with open("{output}.identity_percent_between_BUSCO_duplicated.tsv".format(output=output), 'a') as out:
        out.write("BUSCO_ID\tNumber_of_sequences\tMin_identity\tMax_identity\n")
        for ID, pairs in busco_pairs.items():
            sequences_list, identity_list = [], []
            identity_list.extend([pairs[pair] for pair in pairs])
            for pair in pairs.keys():
                sequences_list.extend(pair.split("_vs_"))

            if len(sequences_list) != 0 and len(identity_list) != 0:
                out.write("{ID}\t{seq}\t{min}\t{max}\n".format(ID=ID, seq=len(set(sequences_list)),
                                                               min=np.min(identity_list), max=np.max(identity_list)))

    with open("{output}.pairs_without_alignments.tsv".format(output=output), 'a') as warning_pairs_output:
        for pair in warning_pairs:
            warning_pairs_output.write("{pair}\n".format(pair=pair))


if __name__ == "__main__":
    blast_dict, busco_dict, busco_pairs, warning_pairs = {}, {}, {}, []
    print("*** BLAST parsing ***")
    blast_parsing(args.blast, blast_dict)
    print("*** BUSCO parsing and pairs creating ***")
    busco_tsv_parsing(args.busco_tsv, busco_dict, busco_pairs, blast_dict, warning_pairs)
    print("*** Output creating ***")
    summary(busco_pairs, warning_pairs, args.output)
