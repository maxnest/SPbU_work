try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--overlaps', type=argparse.FileType('r'), required=True,
                    help="Text file with common GOterms")
parser.add_argument('--externa_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_trunk_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--main_trunk_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--thoracic_part_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_body_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def overlaps_parser(overlaps, overlap_dict):
    for line in overlaps:
        description = line.strip().split("\t")
        samples, terms = description[0], description[1:]
        overlap_dict[samples] = terms

#   print(overlap_dict)


def go_table_parsing(file, go_dict):
    header = file.readline()
    for line in file:
        description = line.strip().split("\t")
        if description[0] not in go_dict.keys():
            go_dict[description[0]] = description[1]


def output_writing(output, overlap_dict, go_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        for samples, terms in overlap_dict.items():
            output_file.write("{samples}\t{annotations}\n".format(
                samples=samples,
                annotations="\t".join(["{term}:{annotation}".format(term=term,
                                                                    annotation=go_dict[term]) for term in terms])))


if __name__ == "__main__":
    overlap_dict, go_dict = {}, {}
    overlaps_parser(args.overlaps, overlap_dict)
    go_table_parsing(args.externa_terms, go_dict)
    go_table_parsing(args.growing_trunk_terms, go_dict)
    go_table_parsing(args.main_trunk_terms, go_dict)
    go_table_parsing(args.thoracic_part_terms, go_dict)
    go_table_parsing(args.whole_body_terms, go_dict)
    output_writing(args.output, overlap_dict, go_dict)