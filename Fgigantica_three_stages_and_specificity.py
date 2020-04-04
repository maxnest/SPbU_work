try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--exp_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--rediae_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--cercariae_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--marita_specific', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def exp_table_parsing(exp_table, gene_dict):
    header = exp_table.readline()
    for line in exp_table:
        description = line.strip().split("\t")
        gene_ID, redia, cercaria, marita = description[0], float(description[1]), \
                                           float(description[2]), float(description[3])

        gene_dict[gene_ID] = {"redia": redia, "cercaria": cercaria, "marita": marita, "specificity": []}


def specificity(gene_dict, table, tag):
    for line in table:
        gene = line.strip()
        gene_dict[gene]["specificity"].append(tag)


def output_writing(output, gene_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tRedia_exp\tCercaria_exp\t"
                          "Marita_exp\tSpecificity\n")
        for gene, values in gene_dict.items():
            output_file.write("{gene_id}\t{redia}\t{cercaria}\t{marita}\t{specificity}\n".format(
                gene_id=gene, redia=values["redia"], cercaria=values["cercaria"], marita=values["marita"],
                specificity=values["specificity"][0]
            ))


if __name__ == "__main__":
    gene_dict = {}
    exp_table_parsing(args.exp_table, gene_dict)
    specificity(gene_dict, args.rediae_specific, "R")
    specificity(gene_dict, args.cercariae_specific, "C")
    specificity(gene_dict, args.marita_specific, "M")
    for gene, values in gene_dict.items():
        if len(values["specificity"]) == 0:
            values["specificity"].append('H')
    output_writing(args.output, gene_dict)