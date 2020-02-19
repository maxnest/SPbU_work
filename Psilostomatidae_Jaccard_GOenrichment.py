try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--red_first_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_first_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_first_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_first_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_first_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_first_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_second_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_second_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_second_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_second_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_second_bp', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_second_mf', type=argparse.FileType('r'), required=True)
parser.add_argument('--first_tag', type=str, required=True, help="For example: Psimillimum")
parser.add_argument('--second_tag', type=str, required=True, help="For example: Spseudoglobulus")
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


def go_enrichment_results_comparison(first_go, second_go, first_tag, second_tag, first_sample_tag, second_sample_tag,
                                     category_tag, Jaccard_dict):
    first_name, second_name = "{first_tag}_{first_sample}".format(first_tag=first_tag, first_sample=first_sample_tag), \
                        "{second_tag}_{second_sample}".format(second_tag=second_tag, second_sample=second_sample_tag)
    intersection = set.intersection(*[set(first_go), set(second_go)])
    Jaccard_dict["{first_name}_vs_{second_name}_{category}".format(first_name=first_name, second_name=second_name,
                                                                      category=category_tag)] = {
        "Jaccard": jaccard_similarity(first_go, second_go),
        "Intersection": intersection,
        "{first_name}_specific".format(first_name=first_name): [term for term in first_go if term not in intersection],
        "{second_name}_specific".format(second_name=second_name): [term for term in second_go if term not in intersection]}


def output_files_creating(Jaccard_dict, first_goenrich_results_dict,
                          second_goenrich_results_dict, first_tag, second_tag, output):
    with open("{output}.tsv".format(output=output), 'a') as Jaccard_output:
        Jaccard_output.write("{first_tag}\{second_tag}\tRediae\tCercariae\tMarita\n".format(first_tag=first_tag,
                                                                                            second_tag=second_tag))
        Jaccard_output.write("Rediae\t{red_vs_red_bp}|{red_vs_red_mf}\t{red_vs_cer_bp}|{red_vs_cer_mf}\t"
                             "{red_vs_mar_bp}|{red_vs_mar_mf}\n".format(
                                red_vs_red_bp=Jaccard_dict["{first}_red_first_vs_{second}_red_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                red_vs_red_mf=Jaccard_dict["{first}_red_first_vs_{second}_red_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                red_vs_cer_bp=Jaccard_dict["{first}_red_first_vs_{second}_cer_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                red_vs_cer_mf=Jaccard_dict["{first}_red_first_vs_{second}_cer_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                red_vs_mar_bp=Jaccard_dict["{first}_red_first_vs_{second}_mar_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                red_vs_mar_mf=Jaccard_dict["{first}_red_first_vs_{second}_mar_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"]))
        Jaccard_output.write("Cercariae\t{cer_vs_red_bp}|{cer_vs_red_mf}\t{cer_vs_cer_bp}|{cer_vs_cer_mf}\t"
                             "{cer_vs_mar_bp}|{cer_vs_mar_mf}\n".format(
                                cer_vs_red_bp=Jaccard_dict["{first}_cer_first_vs_{second}_red_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                cer_vs_red_mf=Jaccard_dict["{first}_cer_first_vs_{second}_red_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                cer_vs_cer_bp=Jaccard_dict["{first}_cer_first_vs_{second}_cer_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                cer_vs_cer_mf=Jaccard_dict["{first}_cer_first_vs_{second}_cer_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                cer_vs_mar_bp=Jaccard_dict["{first}_cer_first_vs_{second}_mar_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                cer_vs_mar_mf=Jaccard_dict["{first}_cer_first_vs_{second}_mar_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"]))
        Jaccard_output.write("Marita\t{mar_vs_red_bp}|{mar_vs_red_mf}\t{mar_vs_cer_bp}|{mar_vs_cer_mf}\t"
                             "{mar_vs_mar_bp}|{mar_vs_mar_mf}\n".format(
                                mar_vs_red_bp=Jaccard_dict["{first}_mar_first_vs_{second}_red_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                mar_vs_red_mf=Jaccard_dict["{first}_mar_first_vs_{second}_red_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                mar_vs_cer_bp=Jaccard_dict["{first}_mar_first_vs_{second}_cer_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                mar_vs_cer_mf=Jaccard_dict["{first}_mar_first_vs_{second}_cer_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                mar_vs_mar_bp=Jaccard_dict["{first}_mar_first_vs_{second}_mar_second_bp".format(
                                    first=first_tag, second=second_tag)]["Jaccard"],
                                mar_vs_mar_mf=Jaccard_dict["{first}_mar_first_vs_{second}_mar_second_mf".format(
                                    first=first_tag, second=second_tag)]["Jaccard"]))

    for pair, values in Jaccard_dict.items():
        first_name, second_name = "{first_tag}_{first_sample}_first".format(first_tag=first_tag,
                                                                            first_sample=pair.split("_")[1]), \
                                  "{second_tag}_{second_sample}_second".format(second_tag=second_tag,
                                                                               second_sample=pair.split("_")[5])
        with open("{first_tag}_{first_sample}_vs_{second_tag}_{second_sample}_{category}.tsv".format(
                first_tag=first_tag, first_sample=pair.split("_")[1],
                second_tag=second_tag, second_sample=pair.split("_")[5],
                category=pair.split("_")[-1]), 'a') \
                as pair_output:
            pair_output.write("Intersection (length={inter_len})\t{first}_specific (length={first_len})\t"
                              "{second}_specific (length={second_len})\n".format(
                                inter_len=len(values["Intersection"]),
                                first="{first_tag}_{first_sample}".format(
                                        first_tag=first_tag, first_sample=pair.split("_")[1]),
                                first_len=len(values["{first_name}_specific".format(first_name=first_name)]),
                                second="{second_tag}_{second_sample}".format(
                                         second_tag=second_tag, second_sample=pair.split("_")[5]),
                                second_len=len(values["{second_name}_specific".format(second_name=second_name)])))
            for intersection_term, first_term, second_term in itertools.zip_longest(values["Intersection"],
                values["{first_name}_specific".format(first_name=first_name)],
                values["{second_name}_specific".format(second_name=second_name)], fillvalue="-"):
                full_go_terms = []
                if intersection_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=intersection_term, description=first_goenrich_results_dict[intersection_term]))
                else:
                    full_go_terms.append(intersection_term)

                if first_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=first_term, description=first_goenrich_results_dict[first_term]))
                else:
                    full_go_terms.append(first_term)

                if second_term != "-":
                    full_go_terms.append("{term}:{description}".format(
                        term=second_term, description=second_goenrich_results_dict[second_term]))
                else:
                    full_go_terms.append(second_term)

                pair_output.write("{terms}\n".format(terms="\t".join(full_go_terms)))

        #    pair_output.write("Intersection (length = {inter_len}):\t{intersection}\n"
        #                      "{first}_specific (length = {first_len}):\t{first_specific}\n"
        #                      "{second}_specific (length = {second_len}):\t{second_specific}\n".format(
        #                       inter_len=len(values["Intersection"]), intersection=";".join(values["Intersection"]),
        #                       first="{first_tag}_{first_sample}".format(first_tag=first_tag,
        #                                                                 first_sample=pair.split("_")[1]),
        #                       first_len=len(values["{first_name}_specific".format(first_name=first_name)]),
        #                       first_specific=";".join(values["{first_name}_specific".format(first_name=first_name)]),
        #                       second="{second_tag}_{second_sample}".format(second_tag=second_tag,
        #                                                                    second_sample=pair.split("_")[5]),
        #                       second_len=len(values["{second_name}_specific".format(second_name=second_name)]),
        #                      second_specific=";".join(values["{second_name}_specific".format(second_name=second_name)]
        #                                               )))


if __name__ == "__main__":
    Jaccard_dict, first_goenrich_results_dict, second_goenrich_results_dict = {}, {}, {}
    red_first_bp, red_first_mf, cer_first_bp, cer_first_mf, mar_first_bp, mar_first_mf, red_second_bp, red_second_mf, \
    cer_second_bp, cer_second_mf, mar_second_bp, mar_second_mf = [], [], [], [], [], [], [], [], [], [], [], []
    print("***** Input files parsing *****")
    go_table_parsing(args.red_first_bp, red_first_bp, first_goenrich_results_dict)
    go_table_parsing(args.red_first_mf, red_first_mf, first_goenrich_results_dict)
    go_table_parsing(args.cer_first_bp, cer_first_bp, first_goenrich_results_dict)
    go_table_parsing(args.cer_first_mf, cer_first_mf, first_goenrich_results_dict)
    go_table_parsing(args.mar_first_bp, mar_first_bp, first_goenrich_results_dict)
    go_table_parsing(args.mar_first_mf, mar_first_mf, first_goenrich_results_dict)
    go_table_parsing(args.red_second_bp, red_second_bp, second_goenrich_results_dict)
    go_table_parsing(args.red_second_mf, red_second_mf, second_goenrich_results_dict)
    go_table_parsing(args.cer_second_bp, cer_second_bp, second_goenrich_results_dict)
    go_table_parsing(args.cer_second_mf, cer_second_mf, second_goenrich_results_dict)
    go_table_parsing(args.mar_second_bp, mar_second_bp, second_goenrich_results_dict)
    go_table_parsing(args.mar_second_mf, mar_second_mf, second_goenrich_results_dict)
    print("***** Comparisons performing *****")
    go_enrichment_results_comparison(red_first_bp, red_second_bp, args.first_tag, args.second_tag,
                                     "red_first", "red_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(red_first_mf, red_second_mf, args.first_tag, args.second_tag,
                                     "red_first", "red_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(red_first_bp, cer_second_bp, args.first_tag, args.second_tag,
                                     "red_first", "cer_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(red_first_mf, cer_second_mf, args.first_tag, args.second_tag,
                                     "red_first", "cer_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(red_first_bp, mar_second_bp, args.first_tag, args.second_tag,
                                     "red_first", "mar_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(red_first_mf, mar_second_mf, args.first_tag, args.second_tag,
                                     "red_first", "mar_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_bp, red_second_bp, args.first_tag, args.second_tag,
                                     "cer_first", "red_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_mf, red_second_mf, args.first_tag, args.second_tag,
                                     "cer_first", "red_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_bp, cer_second_bp, args.first_tag, args.second_tag,
                                     "cer_first", "cer_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_mf, cer_second_mf, args.first_tag, args.second_tag,
                                     "cer_first", "cer_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_bp, mar_second_bp, args.first_tag, args.second_tag,
                                     "cer_first", "mar_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(cer_first_mf, mar_second_mf, args.first_tag, args.second_tag,
                                     "cer_first", "mar_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_bp, red_second_bp, args.first_tag, args.second_tag,
                                     "mar_first", "red_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_mf, red_second_mf, args.first_tag, args.second_tag,
                                     "mar_first", "red_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_bp, cer_second_bp, args.first_tag, args.second_tag,
                                     "mar_first", "cer_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_mf, cer_second_mf, args.first_tag, args.second_tag,
                                     "mar_first", "cer_second", "mf", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_bp, mar_second_bp, args.first_tag, args.second_tag,
                                     "mar_first", "mar_second", "bp", Jaccard_dict)
    go_enrichment_results_comparison(mar_first_mf, mar_second_mf, args.first_tag, args.second_tag,
                                     "mar_first", "mar_second", "mf", Jaccard_dict)
    print("***** Output files creating *****")
    output_files_creating(Jaccard_dict, first_goenrich_results_dict, second_goenrich_results_dict,
                          args.first_tag, args.second_tag, args.output)
