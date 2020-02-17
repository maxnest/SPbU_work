try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with results of sequence annotation with eggNOG-mapper")
parser.add_argument('--tab_with_selected', type=argparse.FileType('r'), required=True,
                    help="Table with sequences with normalized expression levels above threshold. "
                         "For instance, which have expression level higher than 1 scaledTPM")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def eggNOG_mapper_parsing(contig_dict, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            contig_ID, annotation = description[0].split(".p")[0], description[1:]
            if len(annotation) >= 5 and len(annotation[5]) != 0:
                contig_dict[contig_ID] = annotation[5]
    print("{num} of sequences have GO-terms".format(num=len(contig_dict.keys())))


def selected_sequences_table(table, selected):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        contig_ID, scaledTPM = description[0], description[1]
        selected.append(contig_ID)
    print("{num} of sequences are selected".format(num=len(selected)))


def output_creating(contig_dict, selected, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        for contig_ID, terms in contig_dict.items():
            if contig_ID in selected:
                output.write("{gene}\t{terms}\n".format(gene=contig_ID, terms=terms))


if __name__ == "__main__":
    contig_dict, selected = {}, []
    eggNOG_mapper_parsing(contig_dict, args.eggnog)
    selected_sequences_table(args.tab_with_selected, selected)
    output_creating(contig_dict, selected, args.out)
