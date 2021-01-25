try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--clusters_merged', type=argparse.FileType('r'), required=True,
                    help="Merged table with co-expressed genes. "
                         "The last column should contain information about cluster")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratr analysis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def cluster_parsing(clusters_merged, clusters_dict):
    header = clusters_merged.readline()
    for line in clusters_merged:
        description = line.strip().split("\t")
        geneID, trait = description[0], description[1]
        if trait not in clusters_dict.keys():
            clusters_dict[trait] = []
        clusters_dict[trait].append(geneID)


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        geneID, phylostrata = description[0], description[-1][1:-1]
        phylostratr_dict[geneID] = phylostrata


def generalization(clusters_dict, phylostratr_dict, summary_dict):
    phylostrates = [value for value in phylostratr_dict.values()]

    for phylostrata in phylostrates:
        summary_dict[phylostrata] = {cluster: [] for cluster in clusters_dict.keys()}

    for cluster, gene_list in clusters_dict.items():
        for geneID in gene_list:
            summary_dict[phylostratr_dict[geneID]][cluster].append(geneID)


def output_writing(output, summary_dict, clusters_dict):
    clusters = [cluster for cluster in clusters_dict.keys()]

    with open("{output}_co-expressed_genes_in_phylostrates_summary.tsv".format(output=output), 'a') as output_file:
        output_file.write("Phylostrates\t{clusters}\n".format(clusters="\t".join(clusters)))
        for phylostrata, cluster_values in summary_dict.items():
            values = ["{len} ({percent}%)".format(len=len(cluster_values[cluster]),
                                                  percent=round((len(cluster_values[cluster]) /
                                                                 len(clusters_dict[cluster])) * 100, 2))
                      for cluster in clusters]
            output_file.write("{phylo}\t{values}\n".format(phylo=phylostrata, values="\t".join(values)))


if __name__ == "__main__":
    cluster_dict, phylostratr_dict, summary_dict = {}, {}, {}
    cluster_parsing(args.clusters_merged, cluster_dict)
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    generalization(cluster_dict, phylostratr_dict, summary_dict)
    output_writing(args.output, summary_dict, cluster_dict)
