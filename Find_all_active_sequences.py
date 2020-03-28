try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


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
        with open("{out}.all_sequences_active_in_the_{sample}.threshold_{threshold}.txt".format(
                out=out, sample=sample, threshold=threshold), 'a') as output:
            for contig, values in contig_dict.items():
                if values[sample] >= threshold:
                    output.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    header, contig_dict = [], {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, header, contig_dict)
    print("***** Specificity measuring *****")
    specificity(header, contig_dict, args.out, args.threshold)