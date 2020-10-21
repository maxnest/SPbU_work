try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--first_tab', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of GO-enrichment analysis for first dataset")
parser.add_argument('--second_tab', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of GO-enrichment analysis for second dataset")
parser.add_argument('--first_tag', type=str, required=True)
parser.add_argument('--second_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, table_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        goid, goterm, trait = description[0], description[1], description[-1][1:-1]
        if trait not in table_dict.keys():
            table_dict[trait] = []
        table_dict[trait].append("{id}|{term}".format(id=goid, term=goterm))


def all_terms_with_description(dict, terms_with_description):
    for trait, goterm_list in dict.items():
        for goterm in goterm_list:
            if goterm.split("|")[0] not in terms_with_description.keys():
                terms_with_description[goterm.split("|")[0]] = goterm.split("|")[1]


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


def GOenrichment_results_comparison(first_dict, second_dict, Jaccard_dict):
    for first_trait in first_dict.keys():
        first_list = [term.split("|")[0] for term in first_dict[first_trait]]
        for second_trait in second_dict.keys():
            second_list = [term.split("|")[0] for term in second_dict[second_trait]]
            # comparison #
            intersection = set.intersection(*[set(first_list), set(second_list)])
            Jaccard_dict["{first}_vs_{second}".format(first=first_trait, second=second_trait)] = {
                "Jaccard": jaccard_similarity(first_list, second_list), "Intersection": intersection,
                "first_specific": [term for term in first_list if term not in intersection],
                "second_specific": [term for term in second_list if term not in intersection]}


def output_writing(output, first_tag, first_dict, second_tag, second_dict, Jaccard_dict, terms_with_description):
    first_traits, second_traits = [trait for trait in first_dict.keys()], [trait for trait in second_dict.keys()]
    with open("{output}.Jaccard_values.tsv".format(output=output), 'a') as jaccard_output:
        jaccard_output.write("{first_tag}\{second_tag}\t{second_traits}\n".format(
            first_tag=first_tag, second_tag=second_tag, second_traits="\t".join(second_traits)))
        for first_trait in first_traits:
            jaccard_values = []
            for second_trait in second_traits:
                jaccard_values.append(str(Jaccard_dict["{first}_vs_{second}".format(first=first_trait,
                                                                                second=second_trait)]["Jaccard"]))
            jaccard_output.write("{first_trait}\t{jaccard_values}\n".format(first_trait=first_trait,
                                                                            jaccard_values="\t".join(jaccard_values)))

    for pair, values in Jaccard_dict.items():
        first_name, second_name = "{first_tag}_{first_trait}".format(first_tag=first_tag,
                                                                     first_trait=pair.split("_vs_")[0]), \
                                  "{second_tag}_{second_trait}".format(second_tag=second_tag,
                                                                       second_trait=pair.split("_vs_")[1])
        with open("{first_name}_vs_{second_name}_comparison_results.tsv".format(
                    first_name=first_name, second_name=second_name), 'a') as pair_output:
            pair_output.write("Intersection (length={inter_len})\t{first}_specific (length={first_len})\t"
                              "{second}_specific (length={second_len})\n".format(
                                inter_len=len(values["Intersection"]), first=first_name,
                                first_len=len(values["first_specific"]),
                                second=second_name, second_len=len(values["second_specific"])))
            for intersection_term, first_term, second_term in itertools.zip_longest(
                    values["Intersection"], values["first_specific"], values["second_specific"], fillvalue="-"):
                full_go_terms = []
                if intersection_term != "-":
                    full_go_terms.append("{term}|{description}".format(
                        term=intersection_term, description=terms_with_description[intersection_term]))
                else:
                    full_go_terms.append(intersection_term)

                if first_term != "-":
                    full_go_terms.append("{term}|{description}".format(
                        term=first_term, description=terms_with_description[first_term]))
                else:
                    full_go_terms.append(first_term)

                if second_term != "-":
                    full_go_terms.append("{term}|{description}".format(
                        term=second_term, description=terms_with_description[second_term]))
                else:
                    full_go_terms.append(second_term)

                pair_output.write("{terms}\n".format(terms="\t".join(full_go_terms)))


if __name__ == "__main__":
    first_dict, second_dict, Jaccard_dict, terms_with_description = {}, {}, {}, {}
    print("***** Input files parsing *****")
    table_parsing(args.first_tab, first_dict)
    table_parsing(args.second_tab, second_dict)
    all_terms_with_description(first_dict, terms_with_description)
    all_terms_with_description(second_dict, terms_with_description)
    print("***** Data analysis *****")
    GOenrichment_results_comparison(first_dict, second_dict, Jaccard_dict)
    print("***** Output writing *****")
    output_writing(args.output, args.first_tag, first_dict, args.second_tag, second_dict,
                   Jaccard_dict, terms_with_description)