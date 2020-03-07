try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--blast_sp1_vs_sp2', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_sp2_vs_sp1', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp1_tag', type=str, required=True)
parser.add_argument('--sp2_tag', type=str, required=True)
args = parser.parse_args()


def BLAST_parsing(BLAST, sp_dict):
    header = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        ID, hit_name = description[0], description[3]
        if hit_name != "no hits":
            sp_dict[ID] = hit_name


def RBBH_pairs_searching(sp1_dict, sp2_dict, pairs_dict):
    pairs = 0
    for sp1_id, sp1_hit in sp1_dict.items():
        if sp1_hit in sp2_dict.keys() and sp2_dict[sp1_hit] == sp1_id:
            pairs_dict[sp1_id] = sp1_hit
            pairs += 1
    print("***** {pairs} pairs of RBBH pairs are founded *****".format(pairs=pairs))


def output_file_creating(pair_dict, sp1_tag, sp2_tag):
    with open("{sp1_tag}_vs_{sp2_tag}.RBBH_pairs.tsv".format(sp1_tag=sp1_tag, sp2_tag=sp2_tag), 'a') as output:
        output.write("{sp1_tag}\t{sp2_tag}\n".format(sp1_tag=sp1_tag, sp2_tag=sp2_tag))
        for sp1_id, sp2_id in pair_dict.items():
            output.write("{sp1_id}\t{sp2_id}\n".format(sp1_id=sp1_id, sp2_id=sp2_id))


if __name__ == "__main__":
    sp1_dict, sp2_dict, pairs_dict = {}, {}, {}
    BLAST_parsing(args.blast_sp1_vs_sp2, sp1_dict)
    BLAST_parsing(args.blast_sp2_vs_sp1, sp2_dict)
    RBBH_pairs_searching(sp1_dict, sp2_dict, pairs_dict)
    output_file_creating(pairs_dict, args.sp1_tag, args.sp2_tag)