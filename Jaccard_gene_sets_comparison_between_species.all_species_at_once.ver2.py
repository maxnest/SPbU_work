try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--rbbh_merged', type=argparse.FileType('r'), required=True,
                    help="Merged table with RBBH pairs with 5 columns:"
                         "Pair_ID First_species_geneID Second_species_geneID First_species_tag Second_species_tag."
                         "One pair per line should be! "
                         "The table should include information about all analyzed species")
parser.add_argument('--geneset_merged', type=argparse.FileType('r'), required=True,
                    help="Merged table with 2 column: GeneID Sample_description."
                         "One gene per line should be! "
                         "The 'sample description' should consist species_tag used in rbbh_merged!"
                         "The table should include information about all analyzed species")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def RBBH_merged_parsing(rbbh_merged, rbbh_dict):
    header = rbbh_merged.readline()
    for line in rbbh_merged:
        description = line.strip().split("\t")
        pair_ID, first_sp_geneID, second_sp_geneID, first_sp, second_sp = description[0], description[1], \
                                                                          description[2], description[3], description[4]
        pair_key = "{first}|{second}".format(first=first_sp, second=second_sp)

        if pair_key not in rbbh_dict.keys():
            rbbh_dict[pair_key] = {first_sp: {}, second_sp: {}}

        rbbh_dict[pair_key][first_sp][first_sp_geneID] = pair_ID
        rbbh_dict[pair_key][second_sp][second_sp_geneID] = pair_ID


def geneset_merged_parsing(geneset_merged, geneset_dict):
    header = geneset_merged.readline()
    for line in geneset_merged:
        description = line.strip().split("\t")
        geneID, sample = description[0], description[1]
        if sample not in geneset_dict.keys():
            geneset_dict[sample] = []
        geneset_dict[sample].append(geneID)


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


def appending_comparison_result(Jaccard_dict, sample_set, other_sample_set, sample_tag, other_sample_tag):
    if sample_tag not in Jaccard_dict.keys():
        Jaccard_dict[sample_tag] = {}

    Jaccard_dict[sample_tag][other_sample_tag] = {
        "Jaccard": jaccard_similarity(sample_set, other_sample_set),
        "Intersection": set.intersection(*[set(sample_set), set(other_sample_set)])}


def gene_sets_comparison(rbbh_dict, geneset_dict, Jaccard_dict):
    samples_shown = []
    for pair_key, pair_value in rbbh_dict.items():
        first_sp, second_sp = pair_key.split("|")[0], pair_key.split("|")[1]

        first_sp_samples = [sample for sample in geneset_dict.keys() if first_sp in sample]
        second_sp_samples = [sample for sample in geneset_dict.keys() if second_sp in sample]

        if first_sp not in samples_shown:
            print("{sp}: {samples}".format(sp=first_sp, samples=";".join(first_sp_samples)))
            samples_shown.append(first_sp)

        if second_sp not in samples_shown:
            print("{sp}: {samples}".format(sp=second_sp, samples=";".join(second_sp_samples)))
            samples_shown.append(second_sp)

        # Between species using sets of RBBH #
        for first_sp_sample in first_sp_samples:
            first_rbbh_set = [pair_value[first_sp][geneID] for geneID in geneset_dict[first_sp_sample]
                              if geneID in pair_value[first_sp].keys()]
            for second_sp_sample in second_sp_samples:
                second_rbbh_set = [pair_value[second_sp][geneID] for geneID in geneset_dict[second_sp_sample]
                                   if geneID in pair_value[second_sp].keys()]

                appending_comparison_result(Jaccard_dict, first_rbbh_set, second_rbbh_set,
                                            first_sp_sample, second_sp_sample)
                appending_comparison_result(Jaccard_dict, second_rbbh_set, first_rbbh_set,
                                            second_sp_sample, first_sp_sample)

        # Within species using sets of geneIDs #
        for first_sp_sample in first_sp_samples:
            for first_sp_other_sample in first_sp_samples:
                appending_comparison_result(Jaccard_dict,
                                            geneset_dict[first_sp_sample], geneset_dict[first_sp_other_sample],
                                            first_sp_sample, first_sp_other_sample)

        for second_sp_sample in second_sp_samples:
            for second_sp_other_sample in second_sp_samples:
                appending_comparison_result(Jaccard_dict,
                                            geneset_dict[second_sp_sample], geneset_dict[second_sp_other_sample],
                                            second_sp_sample, second_sp_other_sample)


def appending_values_in_output(output_file, samples, Jaccard_dict, key):
    for sample in samples:
        values = []
        for other_sample in samples:
            if key == "Jaccard":
                values.append(str(Jaccard_dict[sample][other_sample][key]))
            elif key == "Intersection":
                values.append(str(len(Jaccard_dict[sample][other_sample][key])))

        output_file.write("{sample}\t{values}\n".format(sample=sample, values="\t".join(values)))


def output_writing(output, Jaccard_dict):
    samples = [sample for sample in Jaccard_dict.keys()]

    with open("{output}.Jaccard_scores.tsv".format(output=output), 'a') as output_scores:
        output_scores.write("Samples\t{samples}\n".format(samples="\t".join(samples)))
        appending_values_in_output(output_scores, samples, Jaccard_dict, "Jaccard")

    with open("{output}.Jaccard_intersection.tsv".format(output=output), 'a') as output_intersections:
        output_intersections.write("Samples\t{samples}\n".format(samples="\t".join(samples)))
        appending_values_in_output(output_intersections, samples, Jaccard_dict, "Intersection")


if __name__ == "__main__":
    Jaccard_dict, rbbh_dict, geneset_dict = {}, {}, {}
    print("***** Input files parsing *****")
    RBBH_merged_parsing(args.rbbh_merged, rbbh_dict)
    geneset_merged_parsing(args.geneset_merged, geneset_dict)
    print("***** GeneSets comparison *****")
    gene_sets_comparison(rbbh_dict, geneset_dict, Jaccard_dict)
    print("***** Output files writing *****")
    output_writing(args.output, Jaccard_dict)