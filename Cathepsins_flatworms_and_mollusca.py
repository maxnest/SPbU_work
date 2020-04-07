try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--blast_flatworms', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_mollusca', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def blast_parsing(seq_dict, blast_tab, tag):
    header = blast_tab.readline()
    for line in blast_tab:
        description = line.strip().split("\t")
        seq_id, hit_name = description[0], description[3]
        if seq_id not in seq_dict.keys():
            seq_dict[seq_id] = {"flatworms_hit": [], "flatworms_identity": [], "flatworms_pct_idn": [],
                                "mollusca_hit": [], "mollusca_identity": [], "mollusca_pct_idn": []}

        if hit_name == "no hits":
            seq_dict[seq_id]["{tag}_hit".format(tag=tag)].append("No hit")
            seq_dict[seq_id]["{tag}_identity".format(tag=tag)].append("No hit")
            seq_dict[seq_id]["{tag}_pct_idn".format(tag=tag)].append("No hit")
        else:
            seq_dict[seq_id]["{tag}_hit".format(tag=tag)].append(hit_name)
            seq_dict[seq_id]["{tag}_identity".format(tag=tag)].append(description[-2])
            seq_dict[seq_id]["{tag}_pct_idn".format(tag=tag)].append(description[-1])


def output_creating(seq_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Sequence_ID\tFlatworms:blast_hit\tFlatworms:blast_identity\tFlatworms:blast_pct_idn\t"
                          "Mollusca:blast_hit\tMollusca:blast_identity\tMollusca:blast_pct_idn\n")
        for seq, values in seq_dict.items():
            output_file.write("{seq}\t{flatworms_hit}\t{flatworms_idn}\t{flatworms_pct}\t"
                              "{mollusca_hit}\t{mollusca_idn}\t{mollusca_pct}\n".format(
                                seq=seq, flatworms_hit=values["flatworms_hit"][0],
                                flatworms_idn=values["flatworms_identity"][0],
                                flatworms_pct=values["flatworms_pct_idn"][0], mollusca_hit=values["mollusca_hit"][0],
                                mollusca_idn=values["mollusca_identity"][0], mollusca_pct=values["mollusca_pct_idn"][0]))


if __name__ == "__main__":
    seq_dict = {}
    blast_parsing(seq_dict, args.blast_flatworms, "flatworms")
    blast_parsing(seq_dict, args.blast_mollusca, "mollusca")
    output_creating(seq_dict, args.output)