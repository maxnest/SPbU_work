try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--summary_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--wanted_cols', type=argparse.FileType('r'), required=True,
                    help="Text file with ID of columns (one per line) which are the script must parsing")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def parsing_wanted(wanted_cols, wanted_dict):
    for line in wanted_cols:
        wanted_dict[line.strip()] = []


def read_summary_table(summary_table, summary_dict, wanted_dict):
    header = summary_table.readline().strip().split("\t")
    for line in summary_table:
        description = line.strip().split("\t")
        gene_ID, trans_ID, protein_ID, values = description[0], description[1], description[2], description[3:]
        summary_dict[gene_ID] = {wanted_key: values[header.index(wanted_key) - 3] for wanted_key in wanted_dict.keys()}


def parsing_table(summary_dict, wanted_dict):
    for gene, values in summary_dict.items():
        for wanted_key, value in values.items():
            if value != "-" and value != "No hit":
                wanted_dict[wanted_key].append(gene)


def output_writing(output, wanted_dict):
    all_annotated = []
    with open("{output}.annotation_summary.tsv".format(output=output), 'a') as output_file:
        output_file.write("Database\tNumber of annotated sequences\n")
        for key, seqs in wanted_dict.items():
            output_file.write("{key}\t{number}\n".format(key=key, number=len(set(seqs))))
            all_annotated.extend(seqs)
        output_file.write("*****\t*****\nTotal number\t{total}\n".format(total=len(set(all_annotated))))


if __name__ == "__main__":
    wanted_dict, summary_dict = {}, {}
    parsing_wanted(args.wanted_cols, wanted_dict)
    read_summary_table(args.summary_table, summary_dict, wanted_dict)
    parsing_table(summary_dict, wanted_dict)
    output_writing(args.output, wanted_dict)

