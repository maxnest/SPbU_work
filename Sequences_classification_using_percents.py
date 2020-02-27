try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--percent', type=int, required=True, help="What percentage of the total gene expression must be "
                                                               "in a sample for this gene to be classified as specific "
                                                               "for the sample? For example, 75 or 80")
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


def classification(percent, sum_exp, gene, sample, samples, out):
    if samples[sample] >= (percent * sum_exp)/100:
        with open("{out}.{sample}_specific_genes.txt".format(out=out, sample=sample), 'a') as output:
            output.write("{gene}\n".format(gene=gene))


def dict_parsing(contig_dict, out, percent):
    for gene, samples in contig_dict.items():
        gene_sum_exp = np.sum([samples[sample] for sample in samples])
        for sample in samples:
            classification(percent, gene_sum_exp, gene, sample, samples, out)


if __name__ == "__main__":
    header, contig_dict = [], {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, header, contig_dict)
    print("***** Specificity measuring *****")
    dict_parsing(contig_dict, args.out, args.percent)