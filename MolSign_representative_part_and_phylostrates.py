try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--merged_repr', type=argparse.FileType('r'), required=True,
                    help="Table with GeneIDs belonging to representative part of molecular signature (first column) "
                         "and tag of samples (second column). "
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


def merged_table_parsing(merged_table, merged_dict):
    header = merged_table.readline()
    for line in merged_table:
        description = line.strip().split("\t")
        geneID, sample = description[0], description[1]
        if sample not in merged_dict.keys():
            merged_dict[sample] = []
        merged_dict[sample].append(geneID)


def phylostratigraphy_content(phylostratr_dict, merged_dict, content_dict):
    phylostrates = set([value for value in phylostratr_dict.values()])
    samples = set([key for key in merged_dict.keys()])

    for phylostrate in phylostrates:
        content_dict[phylostrate] = {sample: [] for sample in samples}

    for sample, gene_set in merged_dict.items():
        for geneID in gene_set:
            content_dict[phylostratr_dict[geneID]][sample].append(geneID)


def output_writing(output, content_dict, merged_dict):
    samples = set([key for key in merged_dict.keys()])

    for sample in samples:
        with open("{output}.{sample}.phylosummary_of_repr_part.tsv".format(
                output=output, sample=sample), 'a') as sample_summary:
            sample_summary.write("Phylostratum\tGeneCount\tPercent\tGeneIDs\n")
            for phylostrate, values in content_dict.items():
                sample_summary.write("{phylo}\t{count}\t{percent}\t{genes}\n".format(
                    phylo=phylostrate, count=len(values[sample]),
                    percent=round((len(values[sample])/len(merged_dict[sample])) * 100, 2),
                    genes="{genes}".format(genes=";".join(values[sample]) if len(values[sample]) != 0 else "-")))

    with open("{output}.all_samples.phylosummary_of_repr_part.gene_count.tsv".format(
            output=output), 'a') as all_samples_summary:
        all_samples_summary.write("Phylostratum\t{samples}\n".format(samples="\t".join(samples)))
        for phylostrate, values in content_dict.items():
            all_samples_summary.write("{phylo}\t{counts}\n".format(
                phylo=phylostrate, counts="\t".join([str(len(values[sample])) for sample in samples])))

    with open("{output}.all_samples.phylosummary_of_repr_part.percent.tsv".format(
            output=output), 'a') as all_samples_summary:
        all_samples_summary.write("Phylostratum\t{samples}\n".format(samples="\t".join(samples)))
        for phylostrate, values in content_dict.items():
            all_samples_summary.write("{phylo}\t{percent}\n".format(
                phylo=phylostrate,
                percent="\t".join([str(round((len(values[sample])/len(merged_dict[sample])) * 100, 2))
                                   for sample in samples])))


if __name__ == "__main__":
    phylostratr_dict, merged_dict, content_dict = {}, {}, {}
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    merged_table_parsing(args.merged_repr, merged_dict)
    phylostratigraphy_content(phylostratr_dict, merged_dict, content_dict)
    output_writing(args.output, content_dict, merged_dict)