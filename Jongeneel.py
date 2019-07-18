try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()

contig_dict = {}
red_list, cer_list, mar_list = [], [], []


def searching(contig_dict, phase_tag, other_phase_1, other_phase_2, phase_list):
    for contig, values in contig_dict.items():
        if values[phase_tag] >= values[other_phase_1] + values[other_phase_2]:
            phase_list.append(contig)


def write_output(phase_list, phase_tag, output):
    with open("{output}_{phase}_associated".format(output=output, phase=phase_tag), 'a') as output_file:
        for el in set(phase_list):
            output_file.write("{el}\n".format(el=el))


head = args.tab.readline()
for line in args.tab:
    description = line.strip().split("\t")
    ID, red, cer, mar = description[0], float(description[15]), float(description[16]), float(description[17])
    #print("{red}, {cer}, {mar}\n".format(red=red, cer=cer, mar=mar))
    contig_dict[ID] = {"red": red, "cer": cer, "mar": mar}


searching(contig_dict, "red", "cer", "mar", red_list)
searching(contig_dict, "cer", "red", "mar", cer_list)
searching(contig_dict, "mar", "red", "cer", mar_list)

write_output(red_list, "rediae", args.out)
write_output(cer_list, "cercaria", args.out)
write_output(mar_list, "marita", args.out)