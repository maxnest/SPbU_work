try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--first_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--second_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--first_tag', type=str, required=True)
parser.add_argument('--second_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, table_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, trait = description[0][1:-1], description[-1][1:-1]
        if trait not in table_dict.keys():
            table_dict[trait] = []
        table_dict[trait].append(gene)


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


def gene_sets_comparison(first_dict, second_dict, Jaccard_dict):
    for first_trait, first_set in first_dict.items():
        for second_trait, second_set in second_dict.items():
            jaccard_key = "{first}_vs_{second}".format(first=first_trait, second=second_trait)
            intersection = set.intersection(*[set(first_set), set(second_set)])
            Jaccard_dict[jaccard_key] = {"Jaccard": jaccard_similarity(first_set, second_set),
                                         "Intersection": intersection,
                                         "first_specific": [gene for gene in first_set if gene not in intersection],
                                         "second_specific": [gene for gene in second_set if gene not in intersection]}


def output_writing(output, first_dict, second_dict, first_tag, second_tag, Jaccard_dict):
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
        with open("{first_name}_vs_{second_name}_gene_sets_comparison_results.tsv".format(
                    first_name=first_name, second_name=second_name), 'a') as pair_output:
            pair_output.write("Intersection (length={inter_len})\t{first}_specific (length={first_len})\t"
                              "{second}_specific (length={second_len})\n".format(
                                inter_len=len(values["Intersection"]),
                                first=first_name, first_len=len(values["first_specific"]),
                                second=second_name, second_len=len(values["second_specific"])))
            for common_gene, first_gene, second_gene in itertools.zip_longest(
                    values["Intersection"], values["first_specific"], values["second_specific"], fillvalue="-"):
                pair_output.write("{common}\t{first_specific}\t{second_specific}\n".format(
                    common=common_gene, first_specific=first_gene, second_specific=second_gene))


if __name__ == "__main__":
    first_dict, second_dict, Jaccard_dict = {}, {}, {}
    print("***** Input files parsing *****")
    table_parsing(args.first_tab, first_dict)
    table_parsing(args.second_tab, second_dict)
    print("***** Gene sets comparison *****")
    gene_sets_comparison(first_dict, second_dict, Jaccard_dict)
    print("***** Output files writing *****")
    output_writing(args.output, first_dict, second_dict, args.first_tag, args.second_tag, Jaccard_dict)
