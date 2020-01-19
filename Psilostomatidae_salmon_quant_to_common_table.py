try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--red_first', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_second', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_first', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_second', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_first', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_second', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def quant_parsing(contig_dict, sample, sample_tag):
    # Name	Length	EffectiveLength	TPM	NumReads
    header = sample.readline()
    for line in sample:
        description = line.strip().split("\t")
        contig_ID, TPM_value = description[0], float(description[3])
        if contig_ID not in contig_dict.keys():
            contig_dict[contig_ID] = {"red_first": 0, "red_second": 0, "cer_first": 0, "cer_second": 0,
                                      "mar_first": 0, "mar_second": 0}
        contig_dict[contig_ID][sample_tag] += TPM_value


def averaging(contig_dict, dict_with_averaged):
    for contig, values in contig_dict.items():
        dict_with_averaged[contig] = {"red": np.mean([values["red_first"], values["red_second"]]),
                                      "cer": np.mean([values["cer_first"], values["cer_second"]]),
                                      "mar": np.mean([values["mar_first"], values["mar_second"]])}


if __name__ == "__main__":
    contig_dict, dict_with_averaged = {}, {}
    print("***** Input files parsing *****")
    quant_parsing(contig_dict, args.red_first, "red_first")
    quant_parsing(contig_dict, args.red_second, "red_second")
    quant_parsing(contig_dict, args.cer_first, "cer_first")
    quant_parsing(contig_dict, args.cer_second, "cer_second")
    quant_parsing(contig_dict, args.mar_first, "mar_first")
    quant_parsing(contig_dict, args.mar_second, "mar_second")
    averaging(contig_dict, dict_with_averaged)
    print("***** Output files writing *****")
    with open("{out}.all_TPM.tsv".format(out=args.out), 'a') as all_TPM:
        all_TPM.write("Contig_ID\tRediae_first\tRediae_second\tCercariae_first\tCercariae_second\t"
                      "Marita_first\tMarita_second\n")
        for contig, values in contig_dict.items():
            all_TPM.write("{contig}\t{red_first}\t{red_second}\t{cer_first}\t{cer_second}\t{mar_first}\t"
                          "{mar_second}\n".format(contig=contig, red_first=values["red_first"],
                                                  red_second=values["red_second"], cer_first=values["cer_first"],
                                                  cer_second=values["cer_second"], mar_first=values["mar_first"],
                                                  mar_second=values["mar_second"]))

    with open("{out}.averaged_TPM.all.tsv".format(out=args.out), 'a') as all_averaged_TPM:
        all_averaged_TPM.write("Contig_ID\tRediae\tCercariae\tMarita\n")
        for contig, values in dict_with_averaged.items():
            all_averaged_TPM.write("{contig}\t{red}\t{cer}\t{mar}\n".format(
                contig=contig, red=values["red"], cer=values["cer"], mar=values["mar"]
            ))

    with open("{out}.averaged_TPM.more_than_1TPM.tsv".format(out=args.out), 'a') as averaged_TPM_without_low:
        averaged_TPM_without_low.write("Contig_ID\tRediae\tCercariae\tMarita\n")
        for contig, values in dict_with_averaged.items():
            if values["red"] > 1 or values["cer"] > 1 or values["mar"] > 1:
                averaged_TPM_without_low.write("{contig}\t{red}\t{cer}\t{mar}\n".format(
                    contig=contig, red=values["red"], cer=values["cer"], mar=values["mar"]
                ))
