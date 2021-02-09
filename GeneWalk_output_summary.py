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
                         "At least 3 columns should be: Species_ID, Reference_protein_ID, Reference_gene_ID")
parser.add_argument('--output', type=str)
args = parser.parse_args()


def hit_map_parsing(hit_map, hit_map_dict):
    header = hit_map.readline()
    for line in hit_map:
        description = line.strip().split("\t")
        query, subject_protein, subject_gene = description[0], description[1], description[2]
        if subject_gene not in hit_map_dict.keys():
            hit_map_dict[subject_gene] = []
        hit_map_dict[subject_gene].append(query)


def genewalk_results_parsing(genewalk_results, genewalk_dict):
    header = genewalk_results.readline()
    for line in genewalk_results:
        description = line.strip().split(",")
        subject_gene, hgnc_symbol, go_name, go_id, go_domain = description[0], description[1], description[3], \
                                                               description[4], description[-14]

        go_key = "{go_id}:{go_name}".format(go_id=go_id, go_name=go_name)
        if go_domain not in genewalk_dict.keys():
            genewalk_dict[go_domain] = {}

        if go_key not in genewalk_dict[go_domain].keys():
            genewalk_dict[go_domain][go_key] = {}

        if subject_gene not in genewalk_dict[go_domain][go_key].keys():
            genewalk_dict[go_domain][go_key][subject_gene] = {"query": [], "hgnc_symbol": hgnc_symbol}

    # print(genewalk_dict.keys())


def go_terms_filling(genewalk_dict, hit_map_dict):
    for subject_gene, query_list in hit_map_dict.items():
        for go_domain in genewalk_dict.keys():
            for go_key, go_values in genewalk_dict[go_domain].items():
                if subject_gene in go_values.keys():
                    go_values[subject_gene]["query"].extend(query_list)


def output_writing(output, genewalk_dict):
    for go_domain in genewalk_dict.keys():
        with open("{output}.{go_domain}.genewalk_output_summary.tsv".format(
                output=output, go_domain="_".join(go_domain.split(" "))), 'a') as output_file:
            output_file.write("GO_term\t#_genes\tSubject:Query_IDs\n")
            for go_key, go_values in genewalk_dict[go_domain].items():
                output_file.write("{term}\t{counts}\t{ids}\n".format(
                    term=go_key, counts=len(go_values.keys()), ids="\t".join(["{id}:{hgnc}:{query}".format(
                        id=subject_gene, hgnc=value["hgnc_symbol"],
                        query=";".join(value["query"])) for subject_gene, value in go_values.items()])))


if __name__ == "__main__":
    hit_map_dict, genewalk_dict = {}, {}
    hit_map_parsing(args.hit_map, hit_map_dict)
    genewalk_results_parsing(args.genewalk_results, genewalk_dict)
    go_terms_filling(genewalk_dict, hit_map_dict)
    output_writing(args.output, genewalk_dict)