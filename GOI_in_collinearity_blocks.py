try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--sp1_gff3', type=argparse.FileType('r'), required=True,
                    help="GFF3-file for first species only with genes of interest")
parser.add_argument('--sp1_tag', type=str, required=True)
parser.add_argument('--sp2_gff3', type=argparse.FileType('r'), required=True,
                    help="GFF3-file for second species only with gene of interest")
parser.add_argument('--sp2_tag', type=str, required=True)
parser.add_argument('--mcscan', type=argparse.FileType('r'), required=True,
                    help="MCScanX output file with collinearity blocks")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def read_gff3(gff3, sp_list):
    for line in gff3:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            if description[-1].startswith("ID="):
                for el in description[-1].split(";"):
                    if "=gene:" in el:
                        gene = el.split(":")[1]
                        if gene not in sp_list:
                            sp_list.append(gene)


def read_mcscan(mcscan, mcscan_dict):
    for line in mcscan:
        if not line.startswith("#"):
            description = line.strip().split(":")
            number, pair = description[0].split('-')[0], description[1].split(" ")[0]
            sp1, sp2 = pair.split("\t")[1], pair.split("\t")[2]
            if "Alignment_{num}".format(num=number) not in mcscan_dict.keys():
                mcscan_dict["Alignment_{num}".format(num=number)] = {"sp1": [], "sp2": []}
            mcscan_dict["Alignment_{num}".format(num=number)]["sp1"].append(sp1)
            mcscan_dict["Alignment_{num}".format(num=number)]["sp2"].append(sp2)


def goi_in_blocks(mcscan_dict, sp1_list, sp2_list, summary_dict):
    for block, values in mcscan_dict.items():
        sp1_intersection = [gene for gene in values["sp1"] if gene in sp1_list]
        sp2_intersection = [gene for gene in values["sp2"] if gene in sp2_list]
        # what percentage of sequences in the block relates to genes on interest (GOI)?
        summary_dict[block] = {"sp1": round(float(len(sp1_intersection)/len(values["sp1"])) * 100),
                               "sp2": round(float(len(sp2_intersection)/len(values["sp2"])) * 100),
                               "sp1_intersection": sp1_intersection,
                               "sp2_intersection": sp2_intersection}


def write_output(summary_dict, out, sp1_tag, sp2_tag):
    with open("{out}.tab".format(out=out), 'a') as output:
        output.write("NB! Values - what percentage of sequences in the block relates to genes on interest\n"
                     "Block\t{sp1}\t{sp2}\t{sp1}_intersection\t{sp2}_intersection\n".format(sp1=sp1_tag, sp2=sp2_tag))
        for block, values in summary_dict.items():
            output.write("{block}\t{sp1}\t{sp2}\t{sp1_inter}\t{sp2_inter}\n".format(block=block,
                           sp1=values["sp1"], sp2=values["sp2"], sp1_inter=";".join(values["sp1_intersection"]),
                           sp2_inter=";".join(values["sp2_intersection"])))


if __name__ == "__main__":
    sp1_list, sp2_list = [], []
    mcscan_dict, summary_dict = {}, {}
    print("***** Parsing of GFF3-files *****")
    read_gff3(args.sp1_gff3, sp1_list)
    read_gff3(args.sp2_gff3, sp2_list)
    print("***** Parsing MCScanX-output *****")
    read_mcscan(args.mcscan, mcscan_dict)
    print("***** Estimation of GOI percent in blocks *****")
    goi_in_blocks(mcscan_dict, sp1_list, sp2_list, summary_dict)
    print("***** Output files creating *****")
    write_output(summary_dict, args.out, args.sp1_tag, args.sp2_tag)
