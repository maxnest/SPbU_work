try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--averaged_exp', type=argparse.FileType('r'), required=True,
                    help="Table with averaged expression values in TPMs")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr")
parser.add_argument('--threshold', type=int, required=True,
                    help="The threshold for gene activity so that it can be considered significantly expressed")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[protein_ID] = mrca_name


def table_parsing(averaged_exp, exp_dict, list_of_common, threshold):
    stages = averaged_exp.readline().strip().split("\t")[1:]
    print("Stages analyzed: {stages}".format(stages=" ".join(stages)))

    for stage in stages:
        exp_dict[stage] = []

    for line in averaged_exp:
        description = line.strip().split("\t")
        seq_ID, exp_in_stages = description[0], [float(el) for el in description[1:]]

        # genes with stage-specific significant expression:
        if sum(el >= threshold for el in exp_in_stages) == 1:
            for i, el in enumerate(exp_in_stages):
                if el >= threshold:
                    exp_dict[stages[i]].append(seq_ID)

        # genes with significant expression in all analyzed stages:
        if all(el >= threshold for el in exp_in_stages):
            list_of_common.append(seq_ID)


def summarization(phylostratr_dict, exp_dict, list_of_common, summary_dict):
    phylostrates = set([value for value in phylostratr_dict.values()])
    stages = [stage for stage in exp_dict.keys()]

    for phylostrate in phylostrates:
        summary_dict[phylostrate] = {stage: [] for stage in stages}
        summary_dict[phylostrate]["common"] = []

    # some percentage of genes may not be included in any of the philostrates,
    # since only long proteins were analyzed
    summary_dict["Without_phylostrates"] = {stage: [] for stage in stages}
    summary_dict["Without_phylostrates"]["common"] = []

    for stage, seq_IDs in exp_dict.items():
        for seq_ID in seq_IDs:
            if seq_ID in phylostratr_dict.keys():
                summary_dict[phylostratr_dict[seq_ID]][stage].append(seq_ID)
            else:
                summary_dict["Without_phylostrates"][stage].append(seq_ID)

    for seq_ID in list_of_common:
        if seq_ID in phylostratr_dict.keys():
            summary_dict[phylostratr_dict[seq_ID]]["common"].append(seq_ID)
        else:
            summary_dict["Without_phylostrates"]["common"].append(seq_ID)


def output_writing(output, threshold, exp_dict, list_of_common, summary_dict):
    stages = [stage for stage in exp_dict.keys()]
    with open("{output}.{threshold}_TPM.tsv".format(output=output, threshold=threshold), 'a') as output_file:
        output_file.write("Phylostrates\t{stages}\tCommon_genes ({common_genes} genes)\n".format(
            stages="\t".join(["{stage} ({genes} genes)".format(stage=stage, genes=len(exp_dict[stage]))
                              for stage in stages]),
            common_genes=len(list_of_common)))

        for phylostrate, seqs_sets in summary_dict.items():
            values = ["{len} ({percent}%)".format(
                len=len(seqs_sets[stage]), percent=round((len(seqs_sets[stage])/len(exp_dict[stage]))*100, 2))
                      for stage in stages]
            values.append("{len} ({percent}%)".format(
                len=len(seqs_sets["common"]), percent=round((len(seqs_sets["common"])/len(list_of_common))*100, 2)))

            output_file.write("{phylo}\t{values}\n".format(phylo=phylostrate, values="\t".join(values)))


if __name__ == "__main__":
    phylostratr_dict, exp_dict, summary_dict, list_of_common = {}, {}, {}, []
    print("***** Input files parsing *****")
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    table_parsing(args.averaged_exp, exp_dict, list_of_common, args.threshold)
    print("***** Data analysis and summarization *****")
    summarization(phylostratr_dict, exp_dict, list_of_common, summary_dict)
    print("***** Output file writing *****")
    output_writing(args.output, args.threshold, exp_dict, list_of_common, summary_dict)