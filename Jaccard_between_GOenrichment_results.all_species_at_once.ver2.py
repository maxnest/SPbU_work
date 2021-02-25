try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--gsea_merged', type=argparse.FileType('r'), required=True,
                    help="Metged table with GeneOntology enrichment analysis (topGO output file). "
                         "The last column should contain sample description: {species_tag}_{sample_tag} "
                         "for example: Fgig_cercaria or Sman_sporocyst")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def gsea_parsing(gsea_merged, gsea_dict):
    header = gsea_merged.readline()
    for line in gsea_merged:
        description = line.strip().split("\t")
        goid, goterm, sample = description[0], description[1], description[-1]
        if sample not in gsea_dict.keys():
            gsea_dict[sample] = []
        gsea_dict[sample].append(goid)


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


def goid_sets_comparison(gsea_dict, Jaccard_dict):
    all_samples = [sample for sample in gsea_dict.keys()]

    for sample in all_samples:
        if sample not in Jaccard_dict.keys():
            Jaccard_dict[sample] = {}

        for other_sample in all_samples:
            Jaccard_dict[sample][other_sample] = {
                "Jaccard": jaccard_similarity(gsea_dict[sample], gsea_dict[other_sample]),
                "Intersection": set.intersection(*[set(gsea_dict[sample]), set(gsea_dict[other_sample])])}


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
    Jaccard_dict, gsea_dict = {}, {}
    print("***** Input files parsing *****")
    gsea_parsing(args.gsea_merged, gsea_dict)
    print("***** Enriched GOterms sets comparison *****")
    goid_sets_comparison(gsea_dict, Jaccard_dict)
    print("***** Output files writing *****")
    output_writing(args.output, Jaccard_dict)