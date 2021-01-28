try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--molsign', type=argparse.FileType('r'), required=True,
                    help="Table with result of sample`s molecular signature analysis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def molsign_summarization(molsign, summary_dict):
    header = molsign.readline()
    for line in molsign:
        description = line.strip().split("\t")
        geneID, phylostrate, exp_pattern = description[0], description[1], description[2]
        if phylostrate not in summary_dict.keys():
            summary_dict[phylostrate] = {"common_exp": [], "diff_exp": []}
        summary_dict[phylostrate][exp_pattern].append(geneID)


def output_writing(output, summary_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Phylostrates\tExp_pattern\tGeneCounts\n")
        for phylostrate, values in summary_dict.items():
            for pattern, geneIDs in values.items():
                output_file.write("{phylo}\t{pattern}\t{count}\n".format(
                    phylo=phylostrate, pattern=pattern, count=len(geneIDs)
                ))


if __name__ == "__main__":
    summary_dict = {}
    molsign_summarization(args.molsign, summary_dict)
    output_writing(args.output, summary_dict)

