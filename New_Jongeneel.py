try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="Table with normalized (mean for all repeats) expression values")
parser.add_argument('--out', type=str)
args = parser.parse_args()


def parsing(tab, contig_dict, samples):
    header = tab.readline()
    for sample in header.strip().split("\t")[1:]:
        samples.append(sample)

    for line in tab:
        description = line.strip().split("\t")
        contig_ID, values = description[0], [float(value) for value in description[1:]]
        contig_dict[contig_ID] = {"{sample}".format(sample=sample): 0 for sample in samples}
        for sample in samples:
            contig_dict[contig_ID]["{sample}".format(sample=sample)] += values[samples.index(sample)]
            

def jongeneel(contig_dict, results_dict, samples):
    for sample in samples:
        results_dict[sample] = []

    for contig, values in contig_dict.items():
        for key_sample in samples:
            other_samples = [sample for sample in samples if sample != key_sample]
            other_sample_values = 0
            for other_sample in other_samples:
                other_sample_values += values[other_sample]
            if values[key_sample] >= other_sample_values:
                results_dict[key_sample].append(contig)


def write_output(results_dict, out):
    for sample, contigs in results_dict.items():
        with open("{out}.{sample}_specific.txt".format(out=out, sample=sample), 'a') as output_file:
            for contig in contigs:
                output_file.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    contig_dict, results_dict, samples = {}, {}, []
    print("***** Input table parsing *****")
    parsing(args.tab, contig_dict, samples)
    print("***** Sample-specific sequences searching *****")
    jongeneel(contig_dict, results_dict, samples)
    print("***** Output files creating *****")
    write_output(results_dict, args.out)