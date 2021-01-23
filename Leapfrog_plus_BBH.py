try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--bbh', type=argparse.FileType('r'), required=True)
parser.add_argument('--leapfrog', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str)
args = parser.parse_args()


def BestBLASTHits_parsing(bbh, hits_dict):
    head = bbh.readline()
    for line in bbh:
        description = line.strip().split("\t")
        seq_ID, hit_name = description[0], description[3]

        if hit_name != "no hits":

            if seq_ID not in hits_dict.keys():
                hits_dict[seq_ID] = {"hit": [], "identity": [], "pct_identity": [],
                                     "e-value": []}

            hits_dict[seq_ID]["hit"].append("{hit_name}|{hit_description}".format(
                hit_name=hit_name, hit_description=description[4]))
            hits_dict[seq_ID]["identity"].append(description[8])
            hits_dict[seq_ID]["pct_identity"].append(description[9])
            hits_dict[seq_ID]["e-value"].append(description[6])


def leapfrog_parsing(leapfrog, hits_dict):
    head = leapfrog.readline()
    for line in leapfrog:
        description = line.strip().split("\t")
        seq_ID, subject_ID, first_evalue, second_evalue = description[0], description[2], \
                                                          description[3], description[4]
        if seq_ID not in hits_dict.keys():
            hits_dict[seq_ID] = {"hit": [], "identity": [], "pct_identity": [],
                                 "e-value": []}

        hits_dict[seq_ID]["hit"].append(subject_ID)
        hits_dict[seq_ID]["identity"].append("-")
        hits_dict[seq_ID]["pct_identity"].append("-")
        hits_dict[seq_ID]["e-value"].append("{first}|{second}".format(first=first_evalue, second=second_evalue))


def output_writing(output, hits_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Query name:\tHit name:\tE-value:\tIdentity:\tPct Idn:\n")
        for seq_ID, values in hits_dict.items():
            output_file.write("{seq_ID}\t{hit}\t{evalue}\t{identity}\t{pct_idn}\n".format(
                seq_ID=seq_ID, hit=values["hit"][0], evalue=values["e-value"][0],
                identity=values["identity"][0], pct_idn=values["pct_identity"][0]))


if __name__ == "__main__":
    hits_dict = {}
    print("***** Input files parsing *****")
    BestBLASTHits_parsing(args.bbh, hits_dict)
    leapfrog_parsing(args.leapfrog, hits_dict)
    print("***** Output file creating *****")
    output_writing(args.output, hits_dict)
