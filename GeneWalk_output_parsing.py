try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--genewalk_results', type=argparse.FileType('r'), required=True,
                    help="The resulted table created by GeneWalk")
parser.add_argument('--hit_map', type=argparse.FileType('r'), required=True,
                    help="The table with hit map between analyzed species and reference (for example: H.sapiens)."
                         "At least 3 columns should be: Species_seq_ID, Reference_protein_ID, Reference_gene_ID")
parser.add_argument('--wanted', type=argparse.FileType('r'), required=True,
                    help="Text file with wanted geneIDs (one per line) for which information is to be extracted")
parser.add_argument('--global_padj_threshold', type=str,
                    help="This is a FDR significance level (for example, 0.1) "
                         "that will be used as a threshold to classify GO annotations as significant "
                         "or not in this particular experimental context.")
parser.add_argument('--output', type=str)
args = parser.parse_args()


def hit_map_parsing(hit_map, hit_map_dict):
    header = hit_map.readline()
    for line in hit_map:
        description = line.strip().split("\t")
        query_gene, subject_gene = description[0], description[1]
        hit_map_dict[query_gene] = subject_gene


def wanted_parsing(wanted, wanted_list):
    for line in wanted:
        wanted_list.append(line.strip().split("\t")[0])


def genewalk_results_parsing(genewalk_results, global_padj_threshold, wanted_list, hit_map_dict, genewalk_dict):
    wanted_subject_genes = \
        set([hit_map_dict[query_gene] for query_gene in wanted_list if query_gene in hit_map_dict.keys()])

    header = genewalk_results.readline()
    for line in genewalk_results:
        description = line.strip().split(",")
        subject_gene, hgnc_symbol, go_name, go_id, go_domain, global_padj = \
            description[0], description[1], description[3], description[4], description[-14], float(description[-11])
        if subject_gene in wanted_subject_genes and global_padj < float(global_padj_threshold):

            go_key = "{go_id}|{go_name}".format(go_id=go_id, go_name=go_name)

            if go_domain not in genewalk_dict.keys():
                genewalk_dict[go_domain] = {}

            if go_key not in genewalk_dict[go_domain].keys():
                genewalk_dict[go_domain][go_key] = []

            genewalk_dict[go_domain][go_key].append(subject_gene)


def output_summary(output, global_padj_threshold, genewalk_dict):
    for go_domain in genewalk_dict.keys():
        with open("{output}.{go_domain}.{padj}.genewalk_summary.tsv".format(output=output,
                                                                   go_domain="_".join(go_domain.split(" ")),
                                                                   padj=global_padj_threshold), 'a') as summary_output:
            summary_output.write("GO_id\tGO_name\tGeneCount\tSubject_gene_IDs\n")
            for go_key, subject_gene_list in genewalk_dict[go_domain].items():
                summary_output.write("{id}\t{name}\t{count}\t{subject_ids}\n".format(
                    id=go_key.split("|")[0], name=go_key.split("|")[1],
                    count=len(set(subject_gene_list)), subject_ids=";".join(set(subject_gene_list))))


if __name__ == "__main__":
    hit_map_dict, genewalk_dict, wanted_list = {}, {}, []
    print("***** Input files parsing *****")
    hit_map_parsing(args.hit_map, hit_map_dict)
    wanted_parsing(args.wanted, wanted_list)
    genewalk_results_parsing(args.genewalk_results, args.global_padj_threshold,
                             wanted_list, hit_map_dict, genewalk_dict)
    print("***** Output files creating *****")
    output_summary(args.output, args.global_padj_threshold, genewalk_dict)