try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--externa_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_body_terms', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def go_table_parsing(file, list, go_dict):
    head = file.readline()
    for line in file:
        description = line.strip().split("\t")
        list.append(description[0])
        if description[0] not in go_dict.keys():
            go_dict[description[0]] = description[1]


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


def go_enrichment_results_comparison(first_go, second_go, first_sample_tag, second_sample_tag, Jaccard_dict):
    intersection = set.intersection(*[set(first_go), set(second_go)])
    Jaccard_dict["{first_tag}_vs_{second_tag}".format(first_tag=first_sample_tag,
                                                      second_tag=second_sample_tag)] = {
        "Jaccard": jaccard_similarity(first_go, second_go),
        "Intersection": intersection,
        "{first_tag}_specific".format(first_tag=first_sample_tag):
            [term for term in first_go if term not in intersection],
        "{second_tag}_specific".format(second_tag=second_sample_tag):
            [term for term in second_go if term not in intersection]}


def output_files_creating(Jaccard_dict, goenrich_results_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as Jaccard_output_bp:
        Jaccard_output_bp.write("Peltogaster_reticulata\tExterna\tGrowing_stolon\tMiddle_stolon\t"
                                "Terminal_stolon\tWhole_body\n")
        Jaccard_output_bp.write("Externa\t{externa_vs_externa}\t{externa_vs_growing}\t{externa_vs_middle}\t"
                                "{externa_vs_terminal}\t{externa_vs_whole}\n".format(
                                 externa_vs_externa=Jaccard_dict["externa_vs_externa"]["Jaccard"],
                                 externa_vs_growing=Jaccard_dict["externa_vs_growing"]["Jaccard"],
                                 externa_vs_middle=Jaccard_dict["externa_vs_middle"]["Jaccard"],
                                 externa_vs_terminal=Jaccard_dict["externa_vs_terminal"]["Jaccard"],
                                 externa_vs_whole=Jaccard_dict["externa_vs_whole"]["Jaccard"]))
        Jaccard_output_bp.write("Growing_stolon\t{growing_vs_externa}\t{growing_vs_growing}\t{growing_vs_middle}\t"
                                "{growing_vs_terminal}\t{growing_vs_whole}\n".format(
                                 growing_vs_externa=Jaccard_dict["externa_vs_growing"]["Jaccard"],
                                 growing_vs_growing=Jaccard_dict["growing_vs_growing"]["Jaccard"],
                                 growing_vs_middle=Jaccard_dict["growing_vs_middle"]["Jaccard"],
                                 growing_vs_terminal=Jaccard_dict["growing_vs_terminal"]["Jaccard"],
                                 growing_vs_whole=Jaccard_dict["growing_vs_whole"]["Jaccard"]))
        Jaccard_output_bp.write("Middle_stolon\t{middle_vs_externa}\t{middle_vs_growing}\t{middle_vs_middle}\t"
                                "{middle_vs_terminal}\t{middle_vs_whole_body}\n".format(
                                 middle_vs_externa=Jaccard_dict["externa_vs_middle"]["Jaccard"],
                                 middle_vs_growing=Jaccard_dict["growing_vs_middle"]["Jaccard"],
                                 middle_vs_middle=Jaccard_dict["middle_vs_middle"]["Jaccard"],
                                 middle_vs_terminal=Jaccard_dict["middle_vs_terminal"]["Jaccard"],
                                 middle_vs_whole_body=Jaccard_dict["middle_vs_whole"]["Jaccard"]))
        Jaccard_output_bp.write("Terminal_stolon\t{terminal_vs_externa}\t{terminal_vs_growing}\t{terminal_vs_middle}\t"
                                "{terminal_vs_terminal}\t{terminal_vs_whole_body}\n".format(
                                 terminal_vs_externa=Jaccard_dict["externa_vs_terminal"]["Jaccard"],
                                 terminal_vs_growing=Jaccard_dict["growing_vs_terminal"]["Jaccard"],
                                 terminal_vs_middle=Jaccard_dict["middle_vs_terminal"]["Jaccard"],
                                 terminal_vs_terminal=Jaccard_dict["terminal_vs_terminal"]["Jaccard"],
                                 terminal_vs_whole_body=Jaccard_dict["terminal_vs_whole"]["Jaccard"]))
        Jaccard_output_bp.write("Whole_body\t{whole_vs_externa}\t{whole_vs_growing}\t{whole_vs_middle}\t"
                                "{whole_vs_terminal}\t{whole_vs_whole}\n".format(
                                 whole_vs_externa=Jaccard_dict["externa_vs_whole"]["Jaccard"],
                                 whole_vs_growing=Jaccard_dict["growing_vs_whole"]["Jaccard"],
                                 whole_vs_middle=Jaccard_dict["middle_vs_whole"]["Jaccard"],
                                 whole_vs_terminal=Jaccard_dict["terminal_vs_whole"]["Jaccard"],
                                 whole_vs_whole=Jaccard_dict["whole_vs_whole"]["Jaccard"]))

    for pair, values in Jaccard_dict.items():
        with open("Preticulata_ref_{pair}.common_and_specific_GOterms.tsv".format(pair=pair), 'a') as pair_output:
            pair_output.write("Intersection (length={inter_len})\t{first}_specific (length={first_len})\t"
                              "{second}_specific (length={second_len})\n".format(
                               inter_len=len(values["Intersection"]), first=pair.split("_")[0],
                               first_len=len(values["{first_sample}_specific".format(
                                   first_sample=pair.split("_")[0])]),
                               second=pair.split("_")[2],
                               second_len=len(values["{second_sample}_specific".format(
                                   second_sample=pair.split("_")[2])])))
            for intersection_term, first_term, second_term in itertools.zip_longest(
                    values["Intersection"],
                    values["{first_sample}_specific".format(first_sample=pair.split("_")[0])],
                    values["{second_sample}_specific".format(second_sample=pair.split("_")[2])], fillvalue="-"):
                full_go_terms = []
                if intersection_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=intersection_term, description=goenrich_results_dict[intersection_term]))
                else:
                    full_go_terms.append(intersection_term)

                if first_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=first_term, description=goenrich_results_dict[first_term]))
                else:
                    full_go_terms.append(first_term)

                if second_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=second_term, description=goenrich_results_dict[second_term]))
                else:
                    full_go_terms.append(second_term)

                pair_output.write("{terms}\n".format(terms="\t".join(full_go_terms)))


if __name__ == "__main__":
    Jaccard_dict, go_enrichment_dict = {}, {}
    externa_terms, growing_terms, middle_terms, terminal_terms, whole_body_terms = [], [], [], [], []
    print("***** Input files parsing *****")
    go_table_parsing(args.externa_terms, externa_terms, go_enrichment_dict)
    go_table_parsing(args.growing_terms, growing_terms, go_enrichment_dict)
    go_table_parsing(args.middle_terms, middle_terms, go_enrichment_dict)
    go_table_parsing(args.terminal_terms, terminal_terms, go_enrichment_dict)
    go_table_parsing(args.whole_body_terms, whole_body_terms, go_enrichment_dict)
    print("***** Comparisons performing *****")
    # Biological processes
    go_enrichment_results_comparison(externa_terms, externa_terms, "externa", "externa", Jaccard_dict)
    go_enrichment_results_comparison(externa_terms, growing_terms, "externa", "growing", Jaccard_dict)
    go_enrichment_results_comparison(externa_terms, middle_terms, "externa", "middle", Jaccard_dict)
    go_enrichment_results_comparison(externa_terms, terminal_terms, "externa", "terminal", Jaccard_dict)
    go_enrichment_results_comparison(externa_terms, whole_body_terms, "externa", "whole", Jaccard_dict)
    go_enrichment_results_comparison(growing_terms, growing_terms, "growing", "growing", Jaccard_dict)
    go_enrichment_results_comparison(growing_terms, middle_terms, "growing", "middle", Jaccard_dict)
    go_enrichment_results_comparison(growing_terms, terminal_terms, "growing", "terminal", Jaccard_dict)
    go_enrichment_results_comparison(growing_terms, whole_body_terms, "growing", "whole", Jaccard_dict)
    go_enrichment_results_comparison(middle_terms, middle_terms, "middle", "middle", Jaccard_dict)
    go_enrichment_results_comparison(middle_terms, terminal_terms, "middle", "terminal", Jaccard_dict)
    go_enrichment_results_comparison(middle_terms, whole_body_terms, "middle", "whole", Jaccard_dict)
    go_enrichment_results_comparison(terminal_terms, terminal_terms, "terminal", "terminal", Jaccard_dict)
    go_enrichment_results_comparison(terminal_terms, whole_body_terms, "terminal", "whole", Jaccard_dict)
    go_enrichment_results_comparison(whole_body_terms, whole_body_terms, "whole", "whole", Jaccard_dict)
    print("***** Output files creating *****")
    output_files_creating(Jaccard_dict, go_enrichment_dict, args.output)
