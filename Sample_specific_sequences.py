try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--threshold', type=int, required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, header, contig_dict):
    header.extend(table.readline().strip().split("\t"))
    print("Samples: {header}".format(header=" ".join(header[1:])))
    for line in table:
        description = line.strip().split("\t")
        contig_ID, values = description[0], description[1:]
        contig_dict[contig_ID] = {sample: 0 for sample in header[1:]}
        for sample in header[1:]:
            contig_dict[contig_ID][sample] += float(values[header[1:].index(sample)])


def specificity(header, contig_dict, out, threshold):
    for sample in header[1:]:
        other_samples = [other_sample for other_sample in header[1:] if other_sample != sample]
        print("Sample: {sample}; Other samples: {other}".format(sample=sample, other=" ".join(other_samples)))
        with open("{out}.{sample}_specific_set.threshold_{threshold}.txt".format(
                out=out, sample=sample, threshold=threshold), 'a') as output:
            for contig, values in contig_dict.items():
                if values[sample] >= threshold and \
                        np.sum([values[other_sample] for other_sample in other_samples]) == 0:
                    output.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    header, contig_dict = [], {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, header, contig_dict)
    print("***** Specificity measuring *****")
    specificity(header, contig_dict, args.out, args.threshold)