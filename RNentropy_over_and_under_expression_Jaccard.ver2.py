try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import collections

parser = argparse.ArgumentParser()
parser.add_argument('--rbbh', type=argparse.FileType('r'), required=True,
                    help="Table with 1-to-1 reciprocal best BLAST hits (genes) between two species")
parser.add_argument('--first_rnentropy', type=argparse.FileType('r'), required=True,
                    help="Table with significant results of analysis performed using RNentropy for first species")
parser.add_argument('--second_rnentropy', type=argparse.FileType('r'), required=True,
                    help="Table with significant results of analysis performed using RNentropy for second species")
parser.add_argument('--first_tag', type=str, required=True)
parser.add_argument('--second_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def rbbh_parsing(rbbh, rbbh_dict_genes, rbbh_dict_pairs):
    header = rbbh.readline()
    for line in rbbh:
        description = line.strip().split("\t")
        rbbh_ID, first_gene, second_gene = description[0], description[1], description[2]
        rbbh_dict_genes[first_gene] = rbbh_ID
        rbbh_dict_genes[second_gene] = rbbh_ID

        rbbh_dict_pairs[rbbh_ID] = {"first_gene": first_gene, "second_gene": second_gene,
                                    "first_overexp": [], "first_underexp": [],
                                    "second_overexp": [], "second_underexp": [],
                                    "over_similar": [], "under_similar": []}


def RNentropy_parsing(rnentropy, dict_with_overexp, dict_with_underexp, rbbh_dict_genes, rbbh_dict_pairs, tag):
    samples_list = []
    samples_list.extend([sample[1:-1] for sample in rnentropy.readline().strip().split("\t")[3:]])
    # print(samples_list)

    for sample in samples_list:
        dict_with_overexp[sample], dict_with_underexp[sample] = [], []

    for line in rnentropy:
        description = line.strip().split("\t")
        gene, samples = description[0][1:-1], description[3:]
        for index, value in enumerate(samples):
            if value == "1":
                dict_with_overexp[samples_list[index]].append(gene)
                if gene in rbbh_dict_genes.keys():
                    rbbh_dict_pairs[rbbh_dict_genes[gene]]["{tag}_overexp".format(tag=tag)].append(samples_list[index])
            elif value == "-1":
                dict_with_underexp[samples_list[index]].append(gene)
                if gene in rbbh_dict_genes.keys():
                    rbbh_dict_pairs[rbbh_dict_genes[gene]]["{tag}_underexp".format(tag=tag)].append(samples_list[index])


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


def lists_comparison(values, tag):
    if len(values["first_{tag}exp".format(tag=tag)]) != 0 and len(values["second_{tag}exp".format(tag=tag)]) != 0:
        first_samples = [sample for sample in values["first_{tag}exp".format(tag=tag)]]
        second_samples = [sample for sample in values["second_{tag}exp".format(tag=tag)]]
        if collections.Counter(first_samples) == collections.Counter(second_samples):
            values["{tag}_similar".format(tag=tag)].append("*")

        first_samples_splitted = [sample.split("_")[1] for sample in values["first_{tag}exp".format(tag=tag)]]
        second_samples_splitted = [sample.split("_")[1] for sample in values["second_{tag}exp".format(tag=tag)]]
        if collections.Counter(first_samples_splitted) == collections.Counter(second_samples_splitted):
            values["{tag}_similar".format(tag=tag)].append("*")


def rbbh_dict_pairs_analysis(rbbh_dict_pairs):
    for rbbh_ID, values in rbbh_dict_pairs.items():
        lists_comparison(values, "over")
        lists_comparison(values, "under")

    keys = ["first_overexp", "second_overexp", "first_underexp", "second_underexp", "over_similar", "under_similar"]
    for rbbh_ID, values in rbbh_dict_pairs.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append("-")


def RNentropy_results_comparison(rbbh_dict_genes, first_dict, second_dict, Jaccard):
    for first_sample, first_gene_set in first_dict.items():
        rbbh_in_first_sample = [rbbh_dict_genes[gene] for gene in first_gene_set if gene in rbbh_dict_genes.keys()]
        for second_sample, second_gene_set in second_dict.items():
            rbbh_in_second_sample = [rbbh_dict_genes[gene] for gene in second_gene_set if gene in rbbh_dict_genes.keys()]
            Jaccard["{first}_vs_{second}".format(first=first_sample, second=second_sample)] = {
                "Jaccard": jaccard_similarity(rbbh_in_first_sample, rbbh_in_second_sample),
                "Intersection": set.intersection(*[set(rbbh_in_first_sample), set(rbbh_in_second_sample)])
            }


def write_jaccard(output, first_tag, second_tag, first_exp, second_exp, Jaccard, tag):
    with open("{output}.Jaccard_between_{tag}-expressed_gene_sets.tsv".format(
            output=output, tag=tag), 'a') as output_file:
        first_samples = [sample for sample in first_exp.keys()]
        second_samples = [sample for sample in second_exp.keys()]
        output_file.write("{first_tag}\{second_tag}\t{second_samples}\n".format(
            first_tag=first_tag, second_tag=second_tag, second_samples="\t".join(second_samples)))
        for first_sample in first_samples:
            jaccard_values = []
            for second_sample in second_samples:
                jaccard_values.append(str(Jaccard["{first}_vs_{second}".format(first=first_sample,
                                                                               second=second_sample)]["Jaccard"]))
            output_file.write("{first}\t{values}\n".format(first=first_sample, values="\t".join(jaccard_values)))

    with open("{output}.Intersection_lengths_between_{tag}-expressed_gene_sets.tsv".format(
            output=output, tag=tag), 'a') as intersection_lengths_output:
        first_samples = [sample for sample in first_exp.keys()]
        second_samples = [sample for sample in second_exp.keys()]
        intersection_lengths_output.write("{first_tag}\{second_tag}\t{second_samples}\n".format(
            first_tag=first_tag, second_tag=second_tag, second_samples="\t".join(second_samples)))
        for first_sample in first_samples:
            intersection_values = []
            for second_sample in second_samples:
                intersection_values.append(str(len(Jaccard["{first}_vs_{second}".format(
                    first=first_sample, second=second_sample)]["Intersection"])))
            intersection_lengths_output.write("{first}\t{values}\n".format(
                first=first_sample, values="\t".join(intersection_values)))

    with open("{output}.Intersection_between_{tag}-expressed_gene_sets.tsv".format(
            output=output, tag=tag), 'a') as intersection_output:
        intersection_output.write("Sample_pair\tIntersection_length\tIntersection\n")
        for pair, values in Jaccard.items():
            if len(values["Intersection"]) != 0:
                intersection_output.write("{pair}\t{length}\t{intersection}\n".format(
                    pair=pair, length=str(len(values["Intersection"])), intersection=";".join(values["Intersection"])))


def write_rbbh_pairs_with_expression(output, first_tag, second_tag, rbbh_dict_pairs):
    with open("{output}.RBBH_pairs_with_expression.tsv".format(output=output), 'a') as output_file:
        output_file.write("RBBH_pair_ID\t{first}_gene\t{second}_gene\tOver-expression_in_{first}\t"
                          "Over-expression_in_{second}\tUnder-expression_in_{first}\tUnder-expression_in_{second}\t"
                          "Over-expression_pattern_is_similar\tUnder-expression_pattern_is_similar\n".format(
                            first=first_tag, second=second_tag))
        for rbbh_ID, values in rbbh_dict_pairs.items():
            output_file.write("{id}\t{first_gene}\t{second_gene}\t{first_overexp}\t{second_overexp}\t"
                              "{first_underexp}\t{second_underexp}\t{over_similar}\t{under_similar}\n".format(
                                id=rbbh_ID, first_gene=values["first_gene"], second_gene=values["second_gene"],
                                first_overexp=";".join(values["first_overexp"]),
                                second_overexp=";".join(values["second_overexp"]),
                                first_underexp=";".join(values["first_underexp"]),
                                second_underexp=";".join(values["second_underexp"]),
                                over_similar=values["over_similar"][0], under_similar=values["under_similar"][0]))


if __name__ == "__main__":
    rbbh_dict_genes, rbbh_dict_pairs = {}, {}
    first_dict_with_overexp, first_dict_with_underexp = {}, {}
    second_dict_with_overexp, second_dict_with_underexp = {}, {}
    Jaccard_overexp_dict, Jaccard_underexp_dict = {}, {}
    print("***** Input files parsing *****")
    rbbh_parsing(args.rbbh, rbbh_dict_genes, rbbh_dict_pairs)
    RNentropy_parsing(args.first_rnentropy, first_dict_with_overexp, first_dict_with_underexp,
                      rbbh_dict_genes, rbbh_dict_pairs, "first")
    RNentropy_parsing(args.second_rnentropy, second_dict_with_overexp, second_dict_with_underexp,
                      rbbh_dict_genes, rbbh_dict_pairs, "second")
    print("***** Data analysis *****")
    RNentropy_results_comparison(rbbh_dict_genes, first_dict_with_overexp, second_dict_with_overexp,
                                 Jaccard_overexp_dict)
    RNentropy_results_comparison(rbbh_dict_genes, first_dict_with_underexp, second_dict_with_underexp,
                                 Jaccard_underexp_dict)
    rbbh_dict_pairs_analysis(rbbh_dict_pairs)
    print("***** Output files creating *****")
    write_jaccard(args.output, args.first_tag, args.second_tag, first_dict_with_overexp, second_dict_with_overexp,
                  Jaccard_overexp_dict, "over")
    write_jaccard(args.output, args.first_tag, args.second_tag, first_dict_with_underexp, second_dict_with_underexp,
                  Jaccard_underexp_dict, "under")
    write_rbbh_pairs_with_expression(args.output, args.first_tag, args.second_tag, rbbh_dict_pairs)