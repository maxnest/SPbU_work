try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--minlen', type=int, required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def outfmt6_parsing(tab, minlen, dict):
    for line in tab:
        description = line.strip().split("\t")
        qseqid, sseqid, qframe, sframe, qstart, qend, sstart, send, pident, length = description[0],\
           description[1], description[2], description[3], description[4], description[5], \
           description[6], description[7], description[8], int(description[9])
        if qseqid != sseqid:
            if length >= minlen:
                if "-" in qframe and "-" not in sframe or "-" not in qframe and "-" in sframe:
                    if "{qseqid}_vs_{sseqid}".format(qseqid=qseqid, sseqid=sseqid) not in dict.keys() and  \
                            "{sseqid}_vs_{qseqid}".format(sseqid=sseqid, qseqid=qseqid) not in dict.keys():
                        dict["{qseqid}_vs_{sseqid}".format(qseqid=qseqid, sseqid=sseqid)] = {"qseqid": qseqid,
                            "sseqid": sseqid, "qframe": qframe, "sframe": sframe, "qstart": qstart, "qend": qend,
                            "sstart": sstart, "send": send, "pident": pident, "length": length}


def write_tab(out, dict):
    with open("{out}.tsv".format(out=out), 'a') as output:
        output.write("QseqID\tSseqID\tQframe\tSframe\tQstart\tQend\tSstart\tSend\tPident\tLength\n")
        for pair, values in dict.items():
          output.write("{qseqid}\t{sseqid}\t{qframe}\t{sframe}\t{qstart}\t{qend}\t{sstart}\t{send}\t"
                       "{pident}\t{length}\n".format(qseqid=values["qseqid"], sseqid=values["sseqid"],
                                                     qframe=values["qframe"], sframe=values["sframe"],
                                                     qstart=values["qstart"], qend=values["qend"],
                                                     sstart=values["sstart"], send=values["send"],
                                                     pident=values["pident"], length=values["length"]))


if __name__ == "__main__":
    blast_tab = {}
    print("***** BLAST parsing *****")
    outfmt6_parsing(args.tab, args.minlen, blast_tab)
    print("***** Output creating ******")
    write_tab(args.out, blast_tab)

