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


def eggNOG_mapper_parsing(list_of_terms, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein_ID, annotation = description[0], description[1:]
            if len(annotation) >= 5 and len(annotation[5]) != 0:
                list_of_terms.extend(annotation[5].split(","))


def output_creating(list_of_terms, out):
    print("In total {length} unique GeneOntology terms were finded".format(length=len(set(list_of_terms))))
    with open("{out}.tsv".format(out=out), 'a') as output:
        for term in set(list_of_terms):
            output.write("{term}\n".format(term=term))


if __name__ == "__main__":
    list_of_terms = []
    eggNOG_mapper_parsing(list_of_terms, args.eggnog)
    output_creating(list_of_terms, args.out)
