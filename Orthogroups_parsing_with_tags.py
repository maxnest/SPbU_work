try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="OrthoFinder output file: Orthogroups.csv")
parser.add_argument('--wanted', type=argparse.FileType('r'), required=True,
                    help="Text file with species names which you want to include in outfiles;"
                         "one species name per line")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def wanted_parsing(wanted, wanted_list):
    for line in wanted:
        wanted_list.append(line.strip())


def orthogroups_parsing(orthogroups, output_dict, wanted_list):
    head = orthogroups.readline().strip().split("\t")
    # print(head)
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteomes = description[0], description[1:]
        if len(proteomes) > max([head.index(wanted_species) for wanted_species in wanted_list]):
            wanted_proteomes = [proteomes[head.index(wanted_species)] for wanted_species in wanted_list]
            # print(wanted_proteomes)
            lengths = 0
            for proteome in wanted_proteomes:
                lengths += len(proteome)

            if lengths != 0:
                output_dict[group_ID] = []
                for proteome in wanted_proteomes:
                    all_sequences = [protein.split(".p")[0] for protein in proteome.split(",")]
                    output_dict[group_ID].append(",".join(all_sequences))


def output_writing(out, wanted_list, output_dict):
    with open("{out}.csv".format(out=out), 'a') as output_file:
        output_file.write("Orthogroups\t{wanted}\n".format(wanted="\t".join(wanted_list)))
        for group_ID, proteomes in output_dict.items():
            output_file.write("{ortho}\t{wanted}\n".format(ortho=group_ID, wanted="\t".join(proteomes)))


if __name__ == "__main__":
    wanted_list, output_dict = [], {}
    wanted_parsing(args.wanted, wanted_list)
    print("Wanted samples: {wanted}".format(wanted=",".join(wanted_list)))
    orthogroups_parsing(args.ortho, output_dict, wanted_list)
    output_writing(args.out, wanted_list, output_dict)

