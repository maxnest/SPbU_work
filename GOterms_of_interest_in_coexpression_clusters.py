try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--gene_map', type=argparse.FileType('r'), required=True,
                    help="Table with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--GOterms', type=argparse.FileType('r'), required=True,
                    help="TSV table with (one per line) "
                         "interesting GO-terms (first column) and description (second column)")
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with annotation of proteins obtained by eggNOG-mapper")
parser.add_argument('--clusters', type=argparse.FileType('r'), required=True,
                    help="Modified table with genes in co-expression clusters")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def gene_map_parsing(gene_map, genes_dict):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_ID, transcript_ID, protein_ID = description[0], description[1], description[2]
        genes_dict[gene_ID] = {"protein_ID": protein_ID, "annotations": []}


def eggNOG_parsing(genes_dict, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein_ID, annotation = description[0], description[1:]
            if len(annotation) >= 5 and len(annotation[5]) != 0:
                for gene_ID, values in genes_dict.items():
                    if protein_ID == values["protein_ID"]:
                        genes_dict[gene_ID]["annotations"].extend(annotation[5].split(","))


def clusters_parsing(clusters, clusters_dict):
    header = clusters.readline().strip().split("\t")

    for el in header:
        clusters_dict[el] = []

    for line in clusters:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                clusters_dict[header[genes.index(gene)]].append(gene)


def GOterms_parsing(GOterms, goterms_dict):
    for line in GOterms:
        description = line.strip().split("\t")
        goterm, term_description = description[0], description[1]
        goterms_dict[goterm] = term_description


def interesting_genes_in_clusters(goterms_dict, genes_dict, clusters_dict, results_dict):
    for term in goterms_dict.keys():
        results_dict[term] = {cluster: [] for cluster in clusters_dict.keys()}
        for gene, values in genes_dict.items():
            if term in values["annotations"]:
                for cluster, coexpressed_genes in clusters_dict.items():
                    if gene in coexpressed_genes:
                        results_dict[term][cluster].append(gene)


def output_writing(output, results_dict, goterms_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        all_terms = [term for term in results_dict.keys()]
        all_clusters = [cluster for cluster in results_dict[all_terms[0]]]
        output_file.write("Clusters\GOterms\t{terms}\n".format(
            terms="\t".join(["{term}|{description}".format(
                term=term, description=goterms_dict[term]) for term in all_terms])))
        for cluster in all_clusters:
            output_file.write("{cluster}\t{numbers}\n".format(
                cluster=cluster,
                numbers="\t".join(["{number}".format(number=len(results_dict[term][cluster])) for term in all_terms])))


if __name__ == "__main__":
    genes_dict, clusters_dict, goterms_dict, results_dict = {}, {}, {}, {}
    gene_map_parsing(args.gene_map, genes_dict)
    eggNOG_parsing(genes_dict, args.eggnog)
    clusters_parsing(args.clusters, clusters_dict)
    GOterms_parsing(args.GOterms, goterms_dict)
    interesting_genes_in_clusters(goterms_dict, genes_dict, clusters_dict, results_dict)
    output_writing(args.output, results_dict, goterms_dict)