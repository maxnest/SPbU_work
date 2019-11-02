try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--first_tab', type=argparse.FileType('r'), required=True,
                    help="Clust output file called Clusters_Objects but without second line")
parser.add_argument('--first_tag', type=str, required=True)
parser.add_argument('--second_tab', type=argparse.FileType('r'), required=True,
                    help="Clust output file called Clusters_Objects but without second line")
parser.add_argument('--second_tag', type=str, required=True)
parser.add_argument('--output', type=str)
args = parser.parse_args()


def read_clusters_tab(table, dict):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)


def jaccard_similarity(one_list, other_list):
    """
    FROM: http://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
    The function measure the similarity between two sets of genes (Jaccard similarity index)
    @param one_list: list of genes
    @param other_list: list of genes
    """
    intersection_cardinality = len(set.intersection(*[set(one_list), set(other_list)]))
    union_cardinality = len(set.union(*[set(one_list), set(other_list)]))
    if union_cardinality != 0:
        return round(intersection_cardinality / float(union_cardinality), 2)
    else:
        return 0


def comparison(first_dict, second_dict, results):
    for cluster, genes in first_dict.items():
        for other_cluster, other_genes in second_dict.items():
            # print("Genes:{genes}\tOther:{other}\n".format(genes="|".join(genes), other="|".join(other_genes)))
            results["{cluster}_vs_{other_cluster}".format(cluster=cluster, other_cluster=other_cluster)] = \
                jaccard_similarity(genes, other_genes)


def output_writing(output, first_dict, first_tag, second_dict, second_tag, results):
    first_dict_keys, second_dict_keys = [key for key in first_dict.keys()], [key for key in second_dict.keys()]
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("{first_tag} clusters\{second_tag} clusters\t{second_clusters}\n".format(
            first_tag=first_tag, second_tag=second_tag, second_clusters="\t".join(second_dict_keys)
        ))
        for cluster in first_dict_keys:
            jaccard_values = []
            for other_cluster in second_dict_keys:
                jaccard_values.append(str(results["{cluster}_vs_{other_cluster}".format(cluster=cluster, other_cluster=other_cluster)]))
            output_file.write("{cluster}\t{jaccard_values}\n".format(cluster=cluster, jaccard_values="\t".join(jaccard_values)))


if __name__ == "__main__":
    first_dict, second_dict, results = {}, {}, {}
    print("***** Input files parsing *****")
    read_clusters_tab(args.first_tab, first_dict)
    read_clusters_tab(args.second_tab, second_dict)
    print("***** Cluster comparison *****")
    comparison(first_dict, second_dict, results)
    print("***** Output file creating *****")
    output_writing(args.output, first_dict, args.first_tag, second_dict, args.second_tag, results)
