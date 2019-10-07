try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def read_tab(samples, tab, dict):
    first_line = tab.readline().strip().split(",")
    for sample in first_line:
        if "IDs" not in sample:
            samples.append("{sample}".format(sample=sample))
    for line in tab:
        description = line.strip().split(",")
        Trans_IDs, counts = description[0], description[1:]
        dict[Trans_IDs] = {"GeneID": Trans_IDs.split(".")[0], "counts": [count for count in counts]}
        #print(dict[Trans_IDs])


def write_output(samples, dict, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        output.write("Transcripts_IDs\tGeneID\t{samples}\n".format(samples="\t".join(samples)))
        for key, values in dict.items():
            output.write("{trans}\t{gene}\t{counts}\n".format(trans=key, gene=values["GeneID"],
                                                              counts="\t".join(values["counts"])))


if __name__ == "__main__":
    tab_dict, samples = {}, []
    print("***** Parsing table *****")
    read_tab(samples, args.tab, tab_dict)
    print("***** Output file creating *****")
    write_output(samples, tab_dict, args.out)
