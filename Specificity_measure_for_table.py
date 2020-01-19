try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, header, contig_dict):
    header.extend(table.readline().strip().split("\t"))
    for line in table:
        description = line.strip().split("\t")
        contig_ID, values = description[0], description[1:]
        contig_dict[contig_ID] = {sample: 0 for sample in header[1:]}
        for sample in header[1:]:
            contig_dict[contig_ID][sample] += float(values[header[1:].index(sample)])


def specificity(header, contig_dict, out):
    for sample in header[1:]:
        other_samples = [other_sample for other_sample in header[1:] if other_sample != sample]
        with open("{out}.{sample}_associated_set.txt".format(out=out, sample=sample), 'a') as output:
            for contig, values in contig_dict.items():
                if values[sample] > np.sum([values[other_sample] for other_sample in other_samples]):
                    output.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    header, contig_dict = [], {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, header, contig_dict)
    print("***** Specificity measuring *****")
    specificity(header, contig_dict, args.out)
