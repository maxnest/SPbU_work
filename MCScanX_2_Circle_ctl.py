try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--collinearity', type=argparse.FileType('r'), required=True)
parser.add_argument('--pixels', type=int, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def synteny_parsing(collinearity, sp1_list, sp2_list):
    for line in collinearity:
        if line.startswith("## Alignment"):
            description = line.strip().split(" ")
            chromosomes = description[6]
            sp1_list.append(chromosomes.split("&")[0])
            sp2_list.append(chromosomes.split("&")[1])


def write_ctl(output, pixels, sp1_list, sp2_list):
    with open("{output}.ctl".format(output=output), 'a') as ctl:
        ctl.write("{pixels}\n".format(pixels=pixels))
        ctl.write("{sp1},{sp2}".format(sp1=",".join(sp1_list), sp2=",".join(sp2_list)))


if __name__ == "__main__":
    sp1_chromosomes, sp2_chromosomes = [], []
    print("*** MCScanX output file parsing ***")
    synteny_parsing(args.collinearity, sp1_chromosomes, sp2_chromosomes)
    print("*** Circle ctl file creating ***")
    write_ctl(args.output, args.pixels, sp1_chromosomes, sp2_chromosomes)