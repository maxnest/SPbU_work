try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--averaged_exp', type=argparse.FileType('r'), required=True,
                    help="Table with averaged between replicates expression")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr."
                         "Gene IDs should be in table!")
parser.add_argument('--threshold', type=str, required=True,
                    help="The threshold value for gene expression. "
                         "Gene will be included in molucular signature of sample "
                         "if expression of this gene equal or more than this value (in TPM) in considered sample")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def averaged_exp_parsing(averaged_exp, samples, exp_dict):
    samples.extend(averaged_exp.readline().strip().split("\t"))
    print("Samples: {samples}".format(samples=" ".join(samples[1:])))
    for line in averaged_exp:
        description = line.strip().split("\t")
        geneID, values = description[0], description[1:]
        exp_dict[geneID] = {sample: 0 for sample in samples[1:]}
        for sample in samples[1:]:
            exp_dict[geneID][sample] += float(values[samples[1:].index(sample)])


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        geneID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[geneID] = "{ps}:{mrca_name}".format(ps=ps, mrca_name=mrca_name)


def molecular_signatures(exp_dict, threshold, molsign_dict):
    for geneID, values in exp_dict.items():
        for sample, exp in values.items():
            if sample not in molsign_dict.keys():
                molsign_dict[sample] = []

            if exp >= int(threshold):
                molsign_dict[sample].append(geneID)

    print("## Molecular signatures of samples ##")
    for sample, molsign in molsign_dict.items():
        print("# {sample}: {molsign} genes #".format(sample=sample, molsign=len(molsign)))


def common_exp_defining(samples, molsign_dict, common_exp):
    common_exp.extend(set.intersection(*[set(molsign_dict[sample]) for sample in samples[1:]]))

    print("## Number of genes with common expression: {common} ##".format(common=len(common_exp)))


def phylostratigraphic_content(phylostratr_dict, molsign_dict, common_exp, content_dict):
    phylostrates = set([value for value in phylostratr_dict.values()])
    samples = set([key for key in molsign_dict.keys()])

    for phylostrate in phylostrates:
        content_dict[phylostrate] = {sample: {"signature": [], "common_exp": [], "DEG": []} for sample in samples}
        content_dict[phylostrate]["common_exp"] = []

    for sample, geneIDs in molsign_dict.items():
        for geneID in geneIDs:
            content_dict[phylostratr_dict[geneID]][sample]["signature"].append(geneID)
            if geneID in common_exp:
                content_dict[phylostratr_dict[geneID]][sample]["common_exp"].append(geneID)
            else:
                content_dict[phylostratr_dict[geneID]][sample]["DEG"].append(geneID)

    for geneID in common_exp:
        content_dict[phylostratr_dict[geneID]]["common_exp"].append(geneID)


def output_writing(output, threshold, phylostratr_dict, molsign_dict, common_exp, content_dict):
    # Preparation #
    samples = set([key for key in molsign_dict.keys()])

    DEG = {}
    for phylostrate, values in content_dict.items():
        for sample, sample_values in values.items():
            if sample != "common_exp":
                if sample not in DEG.keys():
                    DEG[sample] = []
                DEG[sample].extend(sample_values["DEG"])

    # Summary for each sample: GeneID, Phylostrata, Exp_pattern #
    for sample in samples:
        with open("{output}_{sample}_molsignature_summary.{threshold}TPM.tsv".format(
                output=output, sample=sample, threshold=threshold), 'a') as molsign_summary:
            molsign_summary.write("GeneIDs ({molsign_len} genes)\tPhylostrates\tExpression_pattern\n".format(
                molsign_len=len(molsign_dict[sample])))
            for geneID in molsign_dict[sample]:
                molsign_summary.write("{geneID}\t{phylo}\t{exp}\n".format(
                    geneID=geneID, phylo=phylostratr_dict[geneID],
                    exp="{pattern}".format(pattern="common_exp" if geneID in common_exp else "diff_exp")
                ))

    # Phylostratigraphic affiliation for genes with different expression patterns #
    with open("{output}_phylostratigraphic_affiliation_summary.{threshold}TPM.tsv".format(
            output=output, threshold=threshold), 'a') as affiliation_summary:
        affiliation_summary.write("Phylostrates\tCommon_genes\t{DEG}\n".format(
            DEG="\t".join(["{sample}_DEG".format(sample=sample) for sample in samples])))
        # Percents:
        # Common: percents of genes belong to this phylostrate relative to the total number of common expressed genes
        # DEG: percents of genes belong to this phylostrate relative to the total number of DiffExpressed genes
        for phylostrate, values in content_dict.items():
            DEG_values = ["{phylo} ({percent}%)".format(
                phylo=len(content_dict[phylostrate][sample]["DEG"]),
                percent=round((len(content_dict[phylostrate][sample]["DEG"])/len(DEG[sample]))*100, 2))
                for sample in samples]

            affiliation_summary.write("{phylo}\t{common}\t{DEG_values}\n".format(
                phylo=phylostrate,
                common="{count} ({percent}%)".format(
                    count=len(content_dict[phylostrate]["common_exp"]),
                    percent=round((len(content_dict[phylostrate]["common_exp"])/len(common_exp))*100, 2)),
                DEG_values="\t".join(DEG_values)))


if __name__ == "__main__":
    samples, exp_dict, phylostratr_dict, molsign_dict, common_exp, content_dict = [], {}, {}, {}, [], {}
    print("***** Input files parsing *****")
    averaged_exp_parsing(args.averaged_exp, samples, exp_dict)
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    print("***** Analysis *****")
    molecular_signatures(exp_dict, args.threshold, molsign_dict)
    common_exp_defining(samples, molsign_dict, common_exp)
    phylostratigraphic_content(phylostratr_dict, molsign_dict, common_exp, content_dict)
    print("***** Output writing *****")
    output_writing(args.output, args.threshold, phylostratr_dict, molsign_dict, common_exp, content_dict)
