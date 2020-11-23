try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--go_univ', type=argparse.FileType('r'), required=True,
                    help="Table with results of annotation with eggNOG (GeneOntologyUniverse for topGO).")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr")
parser.add_argument('--ordered_phylostrates', type=argparse.FileType('r'), required=True,
                    help="Text file with ordered phylostrates (one per line)")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[protein_ID] = mrca_name


def go_univ_parsing(go_univ, go_univ_dict):
    for line in go_univ:
        description = line.strip().split("\t")
        protein_ID, go_terms = description[0], description[1]
        for go_term in go_terms.split(","):
            if go_term not in go_univ_dict.keys():
                go_univ_dict[go_term] = []
            go_univ_dict[go_term].append(protein_ID)


def ordered_phylostrates_parsing(ordered_phylostrates, list_with_ordered):
    for line in ordered_phylostrates:
        list_with_ordered.append(line.strip()[1:-1])

    print("Ordered phylostrates: {ordered}".format(ordered=",".join(list_with_ordered)))


def summarization(phylostratr_dict, list_with_ordered, go_univ_dict, summary_dict):
    for go_term, protein_IDs in go_univ_dict.items():
        if go_term not in summary_dict.keys():
            summary_dict[go_term] = {phylostrata: [] for phylostrata in list_with_ordered}

        for protein_ID in protein_IDs:
            summary_dict[go_term][phylostratr_dict[protein_ID]].append(protein_ID)


def output_writing(output, list_with_ordered, summary_dict):
    with open("{output}_phylostrates_partial_contribution_in_GOterms_content.protein_IDs.tsv".format(
            output=output), 'a') as output_with_IDs:
        output_with_IDs.write("GO_ID\t{phylostrates}\n".format(phylostrates="\t".join(list_with_ordered)))
        for go_term, phylostrates_values in summary_dict.items():
            values = []
            for phylostrata in list_with_ordered:
                if len(phylostrates_values[phylostrata]) == 0:
                    values.append("-")
                else:
                    values.append(";".join(phylostrates_values[phylostrata]))
            # values = [";".join(phylostrates_values[phylostrata]) for phylostrata in list_with_ordered]
            output_with_IDs.write("{go_term}\t{values}\n".format(go_term=go_term,
                                                                 values="\t".join(values)))

    with open("{output}_phylostrates_partial_contribution_in_GOterms_content.protein_counts.tsv".format(
            output=output), 'a') as output_with_counts:
        output_with_counts.write("GO_ID\t{phylostrates}\n".format(phylostrates="\t".join(list_with_ordered)))
        for go_term, phylostrates_values in summary_dict.items():
            values = [str(len(phylostrates_values[phylostrata])) for phylostrata in list_with_ordered]
            output_with_counts.write("{go_term}\t{values}\n".format(go_term=go_term,
                                                                    values="\t".join(values)))

    with open("{output}_phylostrates_summarized_contribution_in_GOterms_content.protein_counts.tsv".format(
            output=output), 'a') as output_with_sums:
        output_with_sums.write("GO_ID\t{phylostrates}\n".format(phylostrates="\t".join(list_with_ordered)))
        for go_term, phylostrates_values in summary_dict.items():
            sum, values = 0, []
            for phylostrata in list_with_ordered:
                sum += len(phylostrates_values[phylostrata])
                values.append(str(sum))
            output_with_sums.write("{go_term}\t{values}\n".format(go_term=go_term,
                                                                  values="\t".join(values)))


if __name__ == "__main__":
    phylostratr_dict, go_univ_dict, summary_dict = {}, {}, {}
    list_with_ordered = []
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    go_univ_parsing(args.go_univ, go_univ_dict)
    ordered_phylostrates_parsing(args.ordered_phylostrates, list_with_ordered)
    summarization(phylostratr_dict, list_with_ordered, go_univ_dict, summary_dict)
    output_writing(args.output, list_with_ordered, summary_dict)