try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--rbbh_pairs', type=argparse.FileType('r'), required=True)
parser.add_argument('--first_tag', type=str, required=True, help="For example: Psimillimum")
parser.add_argument('--second_tag', type=str, required=True, help="For example: Spseudoglobulus")
parser.add_argument('--first_gene_map', type=argparse.FileType('r'), required=True,
                    help="Table for first species with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--second_gene_map', type=argparse.FileType('r'), required=True,
                    help="Table for second species with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def rbbh_pairs_parsing(rbbh_pairs, first_tag, second_tag, pair_dict):
    header = rbbh_pairs.readline()
    for line in rbbh_pairs:
        description = line.strip().split("\t")
        pair_ID, sp1_protein, sp2_protein = description[0], description[1], description[2]
        pair_dict[pair_ID] = {"{sp1}_protein".format(sp1=first_tag): sp1_protein,
                              "{sp2}_protein".format(sp2=second_tag): sp2_protein,
                              "{sp1}_gene".format(sp1=first_tag): [],
                              "{sp2}_gene".format(sp2=second_tag): []}


def gene_map_parsing(gene_map, sp_tag, pair_dict):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene, transcript, protein = description[0], description[1], description[2]
        for pair, values in pair_dict.items():
            if protein == values["{sp_tag}_protein".format(sp_tag=sp_tag)]:
                values["{sp_tag}_gene".format(sp_tag=sp_tag)].append(gene)


def output_creating(output, first_tag, second_tag, pair_dict):
    with open("{output}.genes.tsv".format(output=output), 'a') as output_file:
        output_file.write("RBBH_pair_ID\t{first}\t{second}\n".format(first=first_tag, second=second_tag))
        for pair, values in pair_dict.items():
            if len(values["{first}_gene".format(first=first_tag)]) != 0 and \
                    len(values["{second}_gene".format(second=second_tag)]) != 0:
                output_file.write("{pair}\t{first}\t{second}\n".format(
                                    pair=pair, first=values["{first}_gene".format(first=first_tag)][0],
                                    second=values["{second}_gene".format(second=second_tag)][0]))


if __name__ == "__main__":
    pair_dict = {}
    rbbh_pairs_parsing(args.rbbh_pairs, args.first_tag, args.second_tag, pair_dict)
    gene_map_parsing(args.first_gene_map, args.first_tag, pair_dict)
    gene_map_parsing(args.second_gene_map, args.second_tag, pair_dict)
    output_creating(args.output, args.first_tag, args.second_tag, pair_dict)
