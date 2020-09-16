try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with results of sequence annotation with eggNOG-mapper")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def eggNOG_mapper_parsing(protein_dict, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein_ID, annotation = description[0], description[1:]
            if len(annotation) >= 5 and len(annotation[5]) != 0:
                protein_dict[protein_ID] = annotation[5]
    print("{num} of sequences have GO-terms".format(num=len(protein_dict.keys())))


def output_creating(protein_dict, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        for protein, terms in protein_dict.items():
            output.write("{protein}\t{terms}\n".format(protein=protein, terms=terms))


if __name__ == "__main__":
    protein_dict = {}
    eggNOG_mapper_parsing(protein_dict, args.eggnog)
    output_creating(protein_dict, args.out)
