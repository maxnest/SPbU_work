try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--table', type=argparse.FileType('r'), required=True)
parser.add_argument('--colname', type=str, required=True, help="Сolumn name with comma-separated GO-annotation results "
                                                               "For instance: EggNOG:GO_terms")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, header, colname, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        for line in table:
            if not line.startswith("###"):
                if line.startswith("Contig_ID"):
                    header.extend(line.strip().split("\t"))
                else:
                    description = line.strip().split("\t")
                    contig_ID, terms = description[0], description[header.index(colname)]
                    if terms != "-":
                        output.write("{gene}\t{terms}\n".format(gene=contig_ID, terms=terms))


if __name__ == "__main__":
    header = []
    table_parsing(args.table, header, args.colname, args.out)
