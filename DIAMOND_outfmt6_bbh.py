try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--outfmt_6', type=argparse.FileType('r'), required=True,
                    help="Table with BLAST results in 6 output format")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def outfmt6_parsing(outfmt_6, outfmt6_dict):
    for line in outfmt_6:
        description = line.strip().split("\t")
        seq_ID, hit_name, bit_score = description[0].strip().split(" ")[0], \
                                      description[1].strip().split(" ")[0], \
                                      description[-1].strip().split(" ")[0]
        if seq_ID not in outfmt6_dict.keys():
            outfmt6_dict[seq_ID] = {}
        hit_key = "{hit}|{score}".format(hit=hit_name, score=bit_score)
        outfmt6_dict[seq_ID][hit_key] = {"bit_score": float(bit_score),
                                         "line": line}


def output_writing(output, outfmt6_dict):
    with open("{output}.bbh.tsv".format(output=output), 'a') as output_file:
        for key, values in outfmt6_dict.items():
            bit_scores = [hit_value["bit_score"] for hit_value in values.values()]

            for hit_key, hit_value in values.items():
                if hit_value["bit_score"] == np.max(bit_scores):
                    output_file.write("{line}".format(line=hit_value["line"]))
                    break


if __name__ == "__main__":
    outfmt6_dict = {}
    print("***** Input table parsing *****")
    outfmt6_parsing(args.outfmt_6, outfmt6_dict)
    print("***** Output file creating *****")
    output_writing(args.output, outfmt6_dict)








