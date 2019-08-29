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

    for key, value in busco_dict.items():
        pairs = ["_vs_".join(pair) for pair in itertools.combinations(value, r=2)]
        busco_pairs[key] = {pair: blast_dict[pair] for pair in pairs if pair.split("_vs_")[0] != pair.split("_vs_")[1] and pair in blast_dict.keys()}
        for pair in pairs:
            if pair not in blast_dict.keys():
                warning_pairs.append(pair)


def summary(busco_pairs, busco_dict, warning_pairs, output):
    less_50, less_60, less_70, less_80, less_90 = [], [], [], [], []
    with open("{output}.tsv".format(output=output), 'a') as out:
        out.write("ID\t#_of_pairs\tmin_identity\tmax_identity\n")
        for key, values in busco_pairs.items():
            pairs_list = []
            identity_list = []
            for pair, identity in values.items():
                pairs_list.append(pair)
                identity_list.append(identity)
            out.write("{ID}\t{pairs}\t{min}\t{max}\n".format(ID=key, pairs=len(pairs_list),
                                                             min=np.min(identity_list), max=np.max(identity_list)))
            if np.min(identity_list) < 50:
                less_50.append(key)

            if 50 < np.min(identity_list) < 60:
                less_60.append(key)

            if 60 < np.min(identity_list) < 70:
                less_70.append(key)

            if 70 < np.min(identity_list) < 80:
                less_80.append(key)

            if 80 < np.min(identity_list) < 90:
                less_90.append(key)

    with open("{output}.summary.tsv".format(output=output), 'a') as summary:
        summary.write("Warning pairs: \n{pairs}\n".format(pairs="\n".join(warning_pairs)))
        summary.write("Percent of groups with minimal identity between pairs less than:\n"
                      "* 50% : {fifty}\n* 50% < and < 60% : {sixty}\n* 60% < and < 70% : {seventy}\n"
                      "* 70% < and < 80% : {eighty}\n* 80% < and < 90% : {ninety}\n".format(
                        fifty=(len(less_50)/len(busco_dict.keys()))*100,
                        sixty=(len(less_60)/len(busco_dict.keys()))*100,
                        seventy=(len(less_70)/len(busco_dict.keys()))*100,
                        eighty=(len(less_80)/len(busco_dict.keys()))*100,
                        ninety=(len(less_90)/len(busco_dict.keys()))*100
        ))

        print("< 50%: {fifty}, < 60%: {sixty}, < 70%: {seventy}, < 80%: {eighty}, < 90%: {ninety} "
              "and len(busco_dict.keys()) = {keys}".format(fifty=len(less_50), sixty=len(less_60), seventy=len(less_70),
                                                           eighty=len(less_80), ninety=len(less_90),
                                                           keys=len(busco_dict.keys())))


if __name__ == "__main__":
    blast_dict, busco_dict, busco_pairs, warning_pairs = {}, {}, {}, []
    print("*** BLAST parsing ***")
    blast_parsing(args.blast, blast_dict)
    print("*** BUSCO parsing and pairs creating ***")
    busco_tsv_parsing(args.busco_tsv, busco_dict, busco_pairs, blast_dict, warning_pairs)
    print("*** Output creating ***")
    summary(busco_pairs, busco_dict, warning_pairs, args.output)






