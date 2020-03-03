try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--gene_map', type=argparse.FileType('r'), required=True,
                    help="Table with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--clust', type=argparse.FileType('r'), required=True,
                    help="Clust output file with cluster IDs (first row) and genes included in these clusters")
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="eggNOG-mapper output file")
parser.add_argument('--kegg', type=str, required=True, help="The ID of the interesting pathway in the KEGG. \n"
                                                            "For instance: ko04310 or ko04350")
parser.add_argument('--out', type=str, required=True, help="Prefix for output files")
args = parser.parse_args()


def gene_map_parsing(gene_map, gene_dict):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_ID, transcript_ID, protein_ID = description[0], description[1], description[2]
        gene_dict[gene_ID] = {"transcript": transcript_ID, "protein": protein_ID, "cluster": [],
                              "pathways": []}


def clusters(gene_dict, table, dict):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)

    for gene in gene_dict.keys():
        for cluster, genes in dict.items():
            if gene in genes:
                gene_dict[gene]["cluster"].append(cluster)

    for gene, values in gene_dict.items():
        if len(values["cluster"]) == 0:
            values["cluster"].append("-")


def eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein, annotation = description[0], description[1:]
            if len(annotation) >= 9 and len(annotation[8]) != 0:
                if protein in prot_2_gene_dict.keys():
                    gene_dict[prot_2_gene_dict[protein]]["pathways"].extend(annotation[8].split(","))

    for gene, values in gene_dict.items():
        if len(values["pathways"]) == 0:
            values["pathways"].append("-")


def output_writing(out, gene_dict, kegg):
    with open("{out}.{kegg}_pathways_in_clusters.tsv".format(out=out, kegg=kegg), 'a') as output:
        output.write("Gene_ID\tCluster\n")
        for gene, values in gene_dict.items():
            if kegg in values["pathways"]:
                output.write("{gene}\t{cluster}\n".format(gene=gene, cluster=values["cluster"][0]))


if __name__ == "__main__":
    gene_dict, cluster_dict = {}, {}
    gene_map_parsing(args.gene_map, gene_dict)
    prot_2_gene_dict = {gene_dict[gene]["protein"]: gene for gene in gene_dict.keys()}
    clusters(gene_dict, args.clust, cluster_dict)
    eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, args.eggnog)
    output_writing(args.out, gene_dict, args.kegg)