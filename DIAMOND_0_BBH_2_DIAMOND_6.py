try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--outfmt_0_bbh', type=argparse.FileType('r'), required=True,
                    help="Table with Best BLAST hits (outfmt 0)")
parser.add_argument('--outfmt_6', type=argparse.FileType('r'), required=True,
                    help="Table with BLAST results in 6 output format")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def bbh_parsing(outfmt_0_bbh, bbh_dict):
    head = outfmt_0_bbh.readline()
    for line in outfmt_0_bbh:
        description = line.strip().split("\t")
        seq_ID, hit_name = description[0].strip().split(" ")[0], description[3].strip().split(" ")[0]
        if hit_name != "no hits":
            bbh_dict[seq_ID] = hit_name


def outfmt6_parsing(outfmt_6, bbh_dict, output_list):
    for line in outfmt_6:
        description = line.strip().split("\t")
        seq_ID, hit_name = description[0].strip().split(" ")[0], description[1].strip().split(" ")[0]
        if seq_ID in bbh_dict.keys() and bbh_dict[seq_ID] == hit_name:
            output_list.append(line)


def output_writing(output, output_list):
    with open("{output}.tab".format(output=output), 'a') as output_file:
        for line in output_list:
            output_file.write(line)


if __name__ == "__main__":
    bbh_dict, output_list = {}, []
    print("***** BBH table parsing *****")
    bbh_parsing(args.outfmt_0_bbh, bbh_dict)
    print("***** Outfmt_6 parsing *****")
    outfmt6_parsing(args.outfmt_6, bbh_dict, output_list)
    print("***** Output file writing *****")
    output_writing(args.output, output_list)