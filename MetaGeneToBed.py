try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--metagene', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str)
args = parser.parse_args()

pred_genes = {}

for line in args.metagene:

    if line.startswith("#"):
        contig_name = []
        contig_name.clear()
        contig_name.append(line.strip().split(" ")[1])
        gc = next(args.metagene)
        module = next(args.metagene)
    else:
        desc = line.strip().split("\t")
        if len(desc) > 1:
            start, end, strand, score = desc[0], desc[1], desc[2], desc[4]
            # заменить num=len() на подсчет количества ключей с названием контига в начале
            pred_genes["{contig}|ORF_{num}".format(contig=contig_name[0], num=len(pred_genes) + 1)] = {
                "start": start, "end": end, "strand": strand, "score": score}

with open("{output}.bed".format(output=args.output), "a") as bed:
    for key, value in pred_genes.items():
        keys = key.split("|")
        bed.write("{contig}\t{start}\t{end}\t{ORF_name}\t{score}\t{strand}\n".format(
            contig=keys[0], start=value["start"], end=value["end"], ORF_name=keys[1],
            score=value["score"], strand=value["strand"]))
