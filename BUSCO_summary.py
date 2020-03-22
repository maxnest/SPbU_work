try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--csinensis', type=argparse.FileType('r'), required=True)
parser.add_argument('--fgigantica', type=argparse.FileType('r'), required=True)
parser.add_argument('--fhepatica', type=argparse.FileType('r'), required=True)
parser.add_argument('--ofelineus', type=argparse.FileType('r'), required=True)
parser.add_argument('--oviverrini', type=argparse.FileType('r'), required=True)
parser.add_argument('--psimillimum', type=argparse.FileType('r'), required=True)
parser.add_argument('--shaematobium', type=argparse.FileType('r'), required=True)
parser.add_argument('--sjaponicum', type=argparse.FileType('r'), required=True)
parser.add_argument('--smansoni', type=argparse.FileType('r'), required=True)
parser.add_argument('--spseudoglobulus', type=argparse.FileType('r'), required=True)
parser.add_argument('--tregenti', type=argparse.FileType('r'), required=True)
parser.add_argument('--tszidati', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def busco_full_tab_parsing(busco_dict, busco_tab, tag):
    for line in busco_tab:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            ID, status = description[0], description[1]
            if ID not in busco_dict.keys():
                busco_dict[ID] = {"csin": [], "fgig": [], "fhep": [], "ofel": [], "oviv": [], "psim": [],
                                  "shae": [], "sjap": [], "sman": [], "spse": [], "treg": [], "tszi": []}

            busco_dict[ID][tag].append(status)


def output_writing(busco_dict, output):
    keys = ["csin", "fgig", "fhep", "ofel", "oviv", "psim", "shae", "sjap", "sman", "spse", "treg", "tszi"]
    with open("{output}.busco_summary.tsv".format(output=output), 'a') as output_file:
        output_file.write("Ortho_ID\tCsinensis_UP000286415\tFgigantica\tFhepatica_UP000230066\t"
                          "Ofelineus_UP000308267\tOviverrini_UP000054324\tPsimillimum\tShaematobium_UP000054474\t"
                          "Sjaponicum_UP000311919\tSmansoni_UP000008854\tSpseudoglobulus\tTregenti\tTszidati\n")
        for ID, values in busco_dict.items():
            output_file.write("{id}\t{all_status}\n".format(id=ID,
                                                            all_status="\t".join([values[key][0] for key in keys])))

    with open("{output}.common_missing.tsv".format(output=output), 'a') as common_missing:
        for ID, values in busco_dict.items():
            all_status = [values[key][0] for key in keys]
            if len(set(all_status)) == 1 and "Missing" in set(all_status):
                common_missing.write("{id}\n".format(id=ID))


if __name__ == "__main__":
    busco_dict = {}
    print("***** Input files parsing *****")
    busco_full_tab_parsing(busco_dict, args.csinensis, "csin")
    busco_full_tab_parsing(busco_dict, args.fgigantica, "fgig")
    busco_full_tab_parsing(busco_dict, args.fhepatica, "fhep")
    busco_full_tab_parsing(busco_dict, args.ofelineus, "ofel")
    busco_full_tab_parsing(busco_dict, args.oviverrini, "oviv")
    busco_full_tab_parsing(busco_dict, args.psimillimum, "psim")
    busco_full_tab_parsing(busco_dict, args.shaematobium, "shae")
    busco_full_tab_parsing(busco_dict, args.sjaponicum, "sjap")
    busco_full_tab_parsing(busco_dict, args.smansoni, "sman")
    busco_full_tab_parsing(busco_dict, args.spseudoglobulus, "spse")
    busco_full_tab_parsing(busco_dict, args.tregenti, "treg")
    busco_full_tab_parsing(busco_dict, args.tszidati, "tszi")
    print("***** Output file creating *****")
    output_writing(busco_dict, args.output)

