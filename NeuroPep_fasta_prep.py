try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()


def table_parsing(tab, pep_dict):
    for line in tab:
        if not line.startswith("ID"):
            description = line.strip().split("\t")
            pep_ID, seq, family = description[0], description[1].replace(" ", "").replace("/", ""), \
                                  description[4].replace(" ", "_").replace("/", "_")
            pep_dict["{Family}|{ID}".format(Family=family, ID=pep_ID)] = seq
    print("Input table includes {count} sequences".format(count=len(pep_dict.keys())))


def output_writing(out, pep_dict):
    with open("{out}.fasta".format(out=out), 'a') as output:
        for ID, seq in pep_dict.items():
            output.write(">{ID}\n{seq}\n".format(ID=ID, seq=seq))


if __name__ == "__main__":
    pep_dict = {}
    table_parsing(args.tab, pep_dict)
    output_writing(args.out, pep_dict)
