try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--exp_tab', type=argparse.FileType('r'), required=True,
                    help="The prepared table with expression values of selected genes")
parser.add_argument('--go_univ', type=argparse.FileType('r'), required=True,
                    help="The GeneOntology (GO) Universe table with 2 columns (one gene per line): "
                         "GeneID and all GOterms! The table should not contain column names")
parser.add_argument('--wanted_go', type=argparse.FileType('r'), required=True,
                    help="The text file with IDs (one per line) of wanted GO terms")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def exp_tab_parsing(exp_tab, annotation_dict):
    header = exp_tab.readline()
    for line in exp_tab:
        annotation_dict[line.strip().split("\t")[0]] = []


def go_univ_parsing(go_univ, go_univ_dict):
    for line in go_univ:
        geneID, goterms = line.strip().split("\t")[0], line.strip().split("\t")[1]
        go_univ_dict[geneID] = goterms.split(",")


def wanted_go_parsing(wanted_go, wanted_go_dict):
    for line in wanted_go:
        go_id, go_term = line.strip().split("\t")[0], line.strip().split("\t")[1]
        wanted_go_dict[go_id] = go_term


def annotation_dict_filling(annotation_dict, go_univ_dict, wanted_go_dict):
    for geneID, values in annotation_dict.items():
        for wanted_go_id in wanted_go_dict:
            if wanted_go_id in go_univ_dict[geneID]:
                values.append(wanted_go_id)


def output_writing(output, wanted_go_dict, annotation_dict):
    wanted_go_list = [go_id for go_id in wanted_go_dict]

    with open("{output}.pheatmap_annotation.tsv".format(output=output), 'a') as output_file:
        output_file.write("GeneIDs\t{wanted}\n".format(wanted="\t".join(["{id}|{term}".format(
            id=go_id, term=wanted_go_dict[go_id]) for go_id in wanted_go_list])))

        for geneID, values in annotation_dict.items():
            output_values = ["yes" if go_id in values else "no" for go_id in wanted_go_list]
            output_file.write("{geneID}\t{output_values}\n".format(geneID=geneID,
                                                                   output_values="\t".join(output_values)))


if __name__ == "__main__":
    annotation_dict, go_univ_dict, wanted_go_dict = {}, {}, {}
    exp_tab_parsing(args.exp_tab, annotation_dict)
    go_univ_parsing(args.go_univ, go_univ_dict)
    wanted_go_parsing(args.wanted_go, wanted_go_dict)
    annotation_dict_filling(annotation_dict, go_univ_dict, wanted_go_dict)
    output_writing(args.output, wanted_go_dict, annotation_dict)