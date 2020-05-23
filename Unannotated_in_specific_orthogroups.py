try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def tab_parsing(tab, tab_dict):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        ortho, specificity, protein, annotations = description[0], description[1], description[2], description[3:]
        if specificity not in tab_dict.keys():
            tab_dict[specificity] = []
        if (annotations.count("-") + annotations.count("No hit")) == len(annotations):
            tab_dict[specificity].append(protein)


def write_output(output, tab_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_table:
        output_table.write("OGs_specificity\tNumber of unannotated\n")
        for specificity, proteins in tab_dict.items():
            output_table.write("{specificity}\t{proteins}\n".format(specificity=specificity, proteins=len(proteins)))


if __name__ == "__main__":
    tab_dict = {}
    tab_parsing(args.tab, tab_dict)
    write_output(args.output, tab_dict)
    