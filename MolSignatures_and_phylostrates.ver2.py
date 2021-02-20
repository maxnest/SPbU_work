try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--common_exp', type=argparse.FileType('r'), required=True,
                    help="Text file with genes (one per line) with common expression")
parser.add_argument('--over_exp', type=argparse.FileType('r'), required=True,
                    help="Table with over-expressed genes (first column) and tag of samples (second column). "
                         "Only one pair (gene and sample) per line should be!")
parser.add_argument('--specific_exp', type=argparse.FileType('r'), required=True,
                    help="Table with specific-expressed genes (first column) and tag of samples (second column). "
                         "Only one pair (gene and sample) per line should be!")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr."
                         "Gene IDs should be in table!")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        geneID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[geneID] = "{ps}:{mrca_name}".format(ps=ps, mrca_name=mrca_name)


def common_exp_parsing(common_exp, common_exp_list):
    for line in common_exp:
        common_exp_list.append(line.strip().split("\t")[0])


def merged_table_parsing(merged_table, merged_dict):
    header = merged_table.readline()
    for line in merged_table:
        description = line.strip().split("\t")
        geneID, sample = description[0], description[1]
        if sample not in merged_dict.keys():
            merged_dict[sample] = []
        merged_dict[sample].append(geneID)


def append_geneID_to_dict(content_dict, phylostratr_dict, sample, gene_set, gene_set_key):
    for geneID in gene_set:
        content_dict[phylostratr_dict[geneID]][sample][gene_set_key].append(geneID)


def phylostratigraphy_content(phylostratr_dict, common_exp_list, overexp_dict, specific_exp_dict, content_dict):
    phylostrates = set([value for value in phylostratr_dict.values()])
    samples = set([key for key in overexp_dict.keys()])

    for phylostrate in phylostrates:
        content_dict[phylostrate] = {sample: {"over_exp_only": [], "specific_exp_only": [],
                                              "over_and_specific_intersection": [],
                                              "over_and_common_intersection": [],
                                              "common_exp_only": []} for sample in samples}

    for sample in samples:
        over_and_specific_intersection = set.intersection(*[set(overexp_dict[sample]), set(specific_exp_dict[sample])])
        over_and_common_intersection = set.intersection(*[set(overexp_dict[sample]), set(common_exp_list)])
        over_exp_only = [geneID for geneID in overexp_dict[sample]
                         if geneID not in over_and_specific_intersection and geneID not in over_and_common_intersection]
        specific_exp_only = \
            [geneID for geneID in specific_exp_dict[sample] if geneID not in over_and_specific_intersection]
        common_exp_only = [geneID for geneID in common_exp_list if geneID not in over_and_common_intersection]

        append_geneID_to_dict(content_dict, phylostratr_dict, sample,
                              over_and_specific_intersection, "over_and_specific_intersection")
        append_geneID_to_dict(content_dict, phylostratr_dict, sample,
                              over_and_common_intersection, "over_and_common_intersection")
        append_geneID_to_dict(content_dict, phylostratr_dict, sample, over_exp_only, "over_exp_only")
        append_geneID_to_dict(content_dict, phylostratr_dict, sample, specific_exp_only, "specific_exp_only")
        append_geneID_to_dict(content_dict, phylostratr_dict, sample, common_exp_only, "common_exp_only")


def output_writing(output, content_dict):
    samples = []
    for phylostrate, values in content_dict.items():
        samples.extend([sample for sample in values.keys() if sample not in samples])

    pattern_dict = {"over_exp_only": "Over-expression",
                    "specific_exp_only": "Sample-specific expression",
                    "over_and_specific_intersection": "Sample-specific over-expression",
                    "over_and_common_intersection": "Over-expression of common gene",
                    "common_exp_only": "Gene with common expression in all analyzed samples"}

    # Summary tables for samples #
    for sample in samples:
        with open("{output}.{sample}.summary_table.tsv".format(output=output, sample=sample), 'a') as sample_summary:
            sample_summary.write("GeneIDs\tPhylostratum\tExpression_pattern\n")
            for phylostrate, values in content_dict.items():
                for gene_set_name, gene_set_list in values[sample].items():
                    for gene_ID in gene_set_list:
                        sample_summary.write("{id}\t{phylostratum}\t{exp_pattern}\n".format(
                            id=gene_ID, phylostratum=phylostrate, exp_pattern=pattern_dict[gene_set_name]))

    # For bar plots #
    for sample in samples:
        with open("{output}.{sample}.summary_for_barplot.tsv".format(output=output, sample=sample), 'a') as bar_summary:
            bar_summary.write("Phylostratum\tExpression_pattern\tGeneCounts\n")
            for phylostrate, values in content_dict.items():
                for gene_set_name, gene_set_list in values[sample].items():
                    bar_summary.write("{phylostratum}\t{exp_pattern}\t{count}\n".format(
                        phylostratum=phylostrate, exp_pattern=pattern_dict[gene_set_name], count=len(gene_set_list)))

    # Summary table for all stages #
    molsing_of_samples = {}
    for sample in samples:
        molsing_of_samples[sample] = []
        for phylostrate, values in content_dict.items():
            for pattern_key in pattern_dict.keys():
                molsing_of_samples[sample].extend(values[sample][pattern_key])

    with open("{output}.all_samples.summary_table.tsv".format(output=output), 'a') as common_summary:
        common_summary.write("Phylostratum\tExpression_pattern\t{samples}\n".format(
            samples="\t".join(
                ["{sample} ({count} genes)".format(sample=sample, count=len(set(molsing_of_samples[sample])))
                 for sample in samples])))

        for phylostrate, values in content_dict.items():
            for pattern_key, pattern_description in pattern_dict.items():
                common_summary.write("{phylostratum}\t{exp_pattern}\t{counts}\n".format(
                    phylostratum=phylostrate, exp_pattern=pattern_description,
                    counts="\t".join(["{count} ({percent}%)".format(
                        count=len(values[sample][pattern_key]),
                        percent=round((len(values[sample][pattern_key])/len(set(molsing_of_samples[sample]))) * 100, 2))
                        for sample in samples])))


if __name__ == "__main__":
    phylostratr_dict, common_exp_list, overexp_dict, specific_exp_dict, content_dict = {}, [], {}, {}, {}
    print("***** Input files parsing *****")
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    common_exp_parsing(args.common_exp, common_exp_list)
    merged_table_parsing(args.over_exp, overexp_dict)
    merged_table_parsing(args.specific_exp, specific_exp_dict)
    print("***** Analysis *****")
    phylostratigraphy_content(phylostratr_dict, common_exp_list, overexp_dict, specific_exp_dict, content_dict)
    print("***** Output writing *****")
    output_writing(args.output, content_dict)



