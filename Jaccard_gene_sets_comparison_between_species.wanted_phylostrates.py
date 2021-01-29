try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--rbbh', type=argparse.FileType('r'), required=True)
parser.add_argument('--first_tab', type=argparse.FileType('r'), required=True,
                    help="The table with averaged expression values for first species")
parser.add_argument('--second_tab', type=argparse.FileType('r'), required=True,
                    help="The table with averaged expression values for second species")
parser.add_argument('--first_phylostratr', type=argparse.FileType('r'), required=True,
                    help="The table with results of phylostratigraphic analysis (phylostratr) "
                         "of the first species gene set. The gene IDs should be in table!")
parser.add_argument('--second_phylostratr', type=argparse.FileType('r'), required=True,
                    help="The table with results of phylostratigraphic analysis (phylostratr) "
                         "of the second species gene set. The gene IDs should be in table!")
parser.add_argument('--wanted_phylostrates', type=argparse.FileType('r'), required=True,
                    help="The text file with wanted phylostrates (one per line) "
                         "which should be included in the analysis")
parser.add_argument('--first_tag', type=str, required=True)
parser.add_argument('--second_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        gene, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[gene] = mrca_name


def wanted_phylostrates_parsing(wanted_phylostrates, wanted_phylostrates_list):
    for line in wanted_phylostrates:
        wanted_phylostrates_list.append(line.strip().split("\t")[0])

    print("*** Genes related to the next phylostrates would be included into analysis: {wanted} ***".format(
        wanted=";".join(wanted_phylostrates_list)))


def RBBH_pairs_parsing(rbbh, rbbh_dict, first_phylostratr, second_phylostratr, wanted_phylostrates_list):
    header = rbbh.readline()
    for line in rbbh:
        description = line.strip().split("\t")
        pair_ID, first_sp, second_sp = description[0], description[1], description[2]
        if first_phylostratr[first_sp] in wanted_phylostrates_list \
                and second_phylostratr[second_sp] in wanted_phylostrates_list:
            rbbh_dict[pair_ID] = {"first_sp": first_sp, "second_sp": second_sp}


def table_parsing(table, phylostratr_dict, wanted_phylostrates_list, table_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, trait = description[0][1:-1], description[-1][1:-1]
        if trait not in table_dict.keys():
            table_dict[trait] = []

        if phylostratr_dict[gene] in wanted_phylostrates_list:
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


def gene_sets_comparison(rbbh_dict, first_dict, second_dict, Jaccard_dict):
    for first_trait, first_set in first_dict.items():
        for second_trait, second_set in second_dict.items():
            jaccard_key = "{first_trait}_vs_{second_trait}".format(first_trait=first_trait, second_trait=second_trait)
            first_rbbh_set = [pair_ID for pair_ID in rbbh_dict.keys() if rbbh_dict[pair_ID]["first_sp"] in first_set]
            second_rbbh_set = [pair_ID for pair_ID in rbbh_dict.keys() if rbbh_dict[pair_ID]["second_sp"] in second_set]
            intersection = set.intersection(*[set(first_rbbh_set), set(second_rbbh_set)])
            Jaccard_dict[jaccard_key] = {
                "Jaccard": jaccard_similarity(first_rbbh_set, second_rbbh_set), "Intersection": intersection,
                "first_specific": [pair_ID for pair_ID in first_rbbh_set if pair_ID not in intersection],
                "second_specific": [pair_ID for pair_ID in second_rbbh_set if pair_ID not in intersection]}


def output_writing(output, rbbh_dict, first_dict, second_dict, first_tag, second_tag, Jaccard_dict):
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
        with open("{first_name}_vs_{second_name}_rbbh_sets_comparison_results.tsv".format(
                first_name=first_name, second_name=second_name), 'a') as pair_output:
            pair_output.write("Intersection (length={inter_len})\t{first}_specific (length={first_len})\t"
                              "{second}_specific (length={second_len})\n".format(
                                inter_len=len(values["Intersection"]),
                                first=first_name, first_len=len(values["first_specific"]),
                                second=second_name, second_len=len(values["second_specific"])))
            for common_rbbh, first_rbbh, second_rbbh in itertools.zip_longest(
                    values["Intersection"], values["first_specific"], values["second_specific"], fillvalue="-"):
                results = []
                if common_rbbh != "-":
                    results.append("{ID}:{first}_and_{second}".format(
                        ID=common_rbbh, first=rbbh_dict[common_rbbh]["first_sp"],
                        second=rbbh_dict[common_rbbh]["second_sp"]))
                else:
                    results.append(common_rbbh)

                if first_rbbh != "-":
                    results.append("{ID}:{first}".format(
                        ID=first_rbbh, first=rbbh_dict[first_rbbh]["first_sp"]))
                else:
                    results.append(first_rbbh)

                if second_rbbh != "-":
                    results.append("{ID}:{second}".format(
                        ID=second_rbbh, second=rbbh_dict[second_rbbh]["second_sp"]))
                else:
                    results.append(second_rbbh)

                pair_output.write("{results}\n".format(results="\t".join(results)))


if __name__ == "__main__":
    first_phylostratr_dict, second_phylostratr_dict, wanted_phylostrates_list = {}, {}, []
    rbbh_dict, first_dict, second_dict, Jaccard_dict = {}, {}, {}, {}
    print("***** Input files parsing *****")
    phylostratr_parsing(args.first_phylostratr, first_phylostratr_dict)
    phylostratr_parsing(args.second_phylostratr, second_phylostratr_dict)
    wanted_phylostrates_parsing(args.wanted_phylostrates, wanted_phylostrates_list)
    RBBH_pairs_parsing(args.rbbh, rbbh_dict, first_phylostratr_dict, second_phylostratr_dict, wanted_phylostrates_list)
    table_parsing(args.first_tab, first_phylostratr_dict, wanted_phylostrates_list, first_dict)
    table_parsing(args.second_tab, second_phylostratr_dict, wanted_phylostrates_list, second_dict)
    print("***** Gene sets comparison *****")
    gene_sets_comparison(rbbh_dict, first_dict, second_dict, Jaccard_dict)
    print("***** Output files writing *****")
    output_writing(args.output, rbbh_dict, first_dict, second_dict, args.first_tag, args.second_tag, Jaccard_dict)
