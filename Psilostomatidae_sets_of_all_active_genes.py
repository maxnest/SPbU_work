try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--table', type=argparse.FileType('r'), required=True)
parser.add_argument('--threshold', type=int, required=True, help="The minimal expression level")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, rediae, cercariae, marita, threshold):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        contig_ID, red_tpm, cer_tpm, mar_tpm = description[0], float(description[1]), float(description[2]), \
                                               float(description[3])
        if red_tpm > threshold:
            rediae.append(contig_ID)

        if cer_tpm > threshold:
            cercariae.append(contig_ID)

        if mar_tpm > threshold:
            marita.append(contig_ID)


def output_writing(out, list, tag):
    with open("{out}.{tag}_all_active_seqs.tsv".format(out=out, tag=tag), 'a') as output:
        for contig in list:
            output.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    rediae, cercariae, marita = [], [], []
    table_parsing(args.table, rediae, cercariae, marita, args.threshold)
    output_writing(args.out, rediae, "rediae")
    output_writing(args.out, cercariae, "cercariae")
    output_writing(args.out, marita, "marita")