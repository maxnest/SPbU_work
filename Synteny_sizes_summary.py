try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--mcscan', type=argparse.FileType('r'), required=True,
                    help="MCScanX output file with collinearity blocks")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def read_mcscan(mcscan, mcscan_dict):
    for line in mcscan:
        if not line.startswith("#"):
            description = line.strip().split(":")
            number, pair = description[0].split('-')[0], description[1].split(" ")[0]
            sp1, sp2 = pair.split("\t")[1], pair.split("\t")[2]
            if "Alignment_{num}".format(num=number) not in mcscan_dict.keys():
                mcscan_dict["Alignment_{num}".format(num=number)] = []
            mcscan_dict["Alignment_{num}".format(num=number)].append(sp1)


def write_output(mcscan_dict, out):
    with open("{out}.tab".format(out=out), 'a') as output:
        output.write("Group\tSize\n")
        for alignment, values in mcscan_dict.items():
            output.write("{alignment}\t{size}\n".format(alignment=alignment, size=len(values)))


if __name__ == "__main__":
    mcscan_dict = {}
    print("***** Parsing MCScanX-output *****")
    read_mcscan(args.mcscan, mcscan_dict)
    print("***** Output file creating *****")
    write_output(mcscan_dict, args.out)
    