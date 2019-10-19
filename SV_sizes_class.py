try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--ref_tag', type=str, required=True)
parser.add_argument('--sp1_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp1_tag', type=str, required=True)
parser.add_argument('--sp2_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp2_tag', type=str, required=True)
parser.add_argument('--sp3_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp3_tag', type=str, required=True)
parser.add_argument('--sp4_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp4_tag', type=str, required=True)
parser.add_argument('--sp5_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp5_tag', type=str, required=True)
parser.add_argument('--sp6_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp6_tag', type=str, required=True)
parser.add_argument('--sp7_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp7_tag', type=str, required=True)
parser.add_argument('--sp8_tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp8_tag', type=str, required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def classification(tab, tag, ref_dict):
    for line in tab:
        description = line.strip().split("\t")
        ID, length = description[0], int(description[4])
        if length <= 1000:
            ref_dict[tag]["1k"].append(ID)
        elif 1000 < length <= 5000:
            ref_dict[tag]["5k"].append(ID)
        elif 5000 < length <= 10000:
            ref_dict[tag]["10k"].append(ID)
        elif 10000 < length <= 25000:
            ref_dict[tag]["25k"].append(ID)
        elif 25000 < length <= 50000:
            ref_dict[tag]["50k"].append(ID)
        elif length > 50000:
            ref_dict[tag][">50k"].append(ID)


def write_output(ref_dict, out, ref):
    with open("{out}.tab".format(out=out), 'a') as output:
        output.write("{ref}_vs\t0<AND<=1k\t1k<AND<=5k\t5k<AND<=10k\t10k<AND<=25k\t25k<AND<=50k\t>50k\n".format(ref=ref))
        for query, values in ref_dict.items():
            output.write("{query}\t{one_k}\t{five_k}\t{ten_k}\t{twentyfive_k}\t{fifty_k}\t{more_k}\n".format(query=query,
                          one_k=len(values["1k"]), five_k=len(values["5k"]), ten_k=len(values["10k"]),
                          twentyfive_k=len(values["25k"]), fifty_k=len(values["50k"]), more_k=len(values[">50k"])))


if __name__ == "__main__":
    ref_dict = {}
    queries_list = [args.sp1_tag, args.sp2_tag, args.sp3_tag, args.sp4_tag, args.sp5_tag,
                    args.sp6_tag, args.sp7_tag, args.sp8_tag]
    for query in queries_list:
        ref_dict[query] = {"1k": [], "5k": [], "10k": [], "25k": [], "50k": [], ">50k": []}
    print("***** Classification *****")
    classification(args.sp1_tab, args.sp1_tag, ref_dict)
    classification(args.sp2_tab, args.sp2_tag, ref_dict)
    classification(args.sp3_tab, args.sp3_tag, ref_dict)
    classification(args.sp4_tab, args.sp4_tag, ref_dict)
    classification(args.sp5_tab, args.sp5_tag, ref_dict)
    classification(args.sp6_tab, args.sp6_tag, ref_dict)
    classification(args.sp7_tab, args.sp7_tag, ref_dict)
    classification(args.sp8_tab, args.sp8_tag, ref_dict)
    print("***** Output file creating *****")
    write_output(ref_dict, args.out, args.ref_tag)