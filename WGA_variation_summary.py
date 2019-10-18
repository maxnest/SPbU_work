try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--species', type=argparse.FileType('r'), required=True,
                    help="File with species name one per line")
parser.add_argument('--GAP', type=argparse.FileType('r'), required=True,
                    help="File with counts of GAP: fgrep -c 'GAP' */*.show_all_diff.tab")
parser.add_argument('--DUP', type=argparse.FileType('r'), required=True,
                    help="File with counts of DUP: fgrep -c 'DUP' */*.show_all_diff.tab")
parser.add_argument('--BRK', type=argparse.FileType('r'), required=True,
                    help="File with counts of BRK: fgrep -c 'BRK' */*.show_all_diff.tab")
parser.add_argument('--JMP', type=argparse.FileType('r'), required=True,
                    help="File with counts of JMP: fgrep -c 'JMP' */*.show_all_diff.tab")
parser.add_argument('--INV', type=argparse.FileType('r'), required=True,
                    help="File with counts of INV: fgrep -c 'INV' */*.show_all_diff.tab")
parser.add_argument('--SEQ', type=argparse.FileType('r'), required=True,
                    help="File with counts of SEQ: fgrep -c 'SEQ' */*.show_all_diff.tab")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def species_parsing(file, species_list, wga_dict):
    for line in file:
        species_list.append(line.strip())

    for reference in species_list:
        wga_dict[reference] = {"GAP": {}, "DUP": {}, "BRK": {}, "JMP": {}, "INV": {}, "SEQ": {}}
        for query in species_list:
            if query != reference:
                wga_dict[reference]["GAP"][query] = 0
                wga_dict[reference]["DUP"][query] = 0
                wga_dict[reference]["BRK"][query] = 0
                wga_dict[reference]["JMP"][query] = 0
                wga_dict[reference]["INV"][query] = 0
                wga_dict[reference]["SEQ"][query] = 0


def var_parsing(file, var_tag, wga_dict, species_list):
    for line in file:
        description = line.strip().split("/")
        sp_pair, counts = description[0], int(description[-1].split(":")[1])
        ref, query = sp_pair.split("_vs_")[0], sp_pair.split("_vs_")[1]
        if ref in species_list and query in species_list:
            wga_dict[ref][var_tag][query] += counts
        else:
            print("Warning! {ref} or {query} is not in file with species name".format(ref=ref, query=query))


def write_summary(wga_dict, species_list, out):
    with open("{out}.summary.tab".format(out=out), 'a') as summary:
        summary.write("Ref\Queries\t{query}\n".format(query="\t".join(species_list)))
        for reference in species_list:
            results = []
            for query in species_list:
                if reference == query:
                    results.append("X")
                else:
                    results.append("GAP:{GAP}|DUP:{DUP}|BRK:{BRK}|JMP:{JMP}|INV:{INV}|SEQ:{SEQ}".format(
                        GAP=wga_dict[reference]["GAP"][query], DUP=wga_dict[reference]["DUP"][query],
                        BRK=wga_dict[reference]["BRK"][query], JMP=wga_dict[reference]["JMP"][query],
                        INV=wga_dict[reference]["INV"][query], SEQ=wga_dict[reference]["SEQ"][query]
                    ))
            summary.write("{ref}\t{results}\n".format(ref=reference, results="\t".join(results)))

    variance = ["GAP", "DUP", "BRK", "JMP", "INV", "SEQ"]
    for var in variance:
        with open("{out}.{var}.tab".format(out=out, var=var), 'a') as var_output:
            var_output.write("Ref\Quiries\t{query}\n".format(query="\t".join(species_list)))
            for reference in species_list:
                results = []
                for query in species_list:
                    if reference == query:
                        results.append("X")
                    else:
                        results.append(str(wga_dict[reference][var][query]))
                var_output.write("{ref}\t{results}\n".format(ref=reference, results="\t".join(results)))

    with open("{out}.dominants.tab".format(out=out), 'a') as dominant_output:
        dominant_output.write("Ref\Quiries\t{query}\n".format(query="\t".join(species_list)))
        for reference in species_list:
            dominants = []
            for query in species_list:
                if reference == query:
                    dominants.append("X")
                else:
                    var_counts = [wga_dict[reference]["GAP"][query], wga_dict[reference]["DUP"][query],
                                  wga_dict[reference]["BRK"][query], wga_dict[reference]["JMP"][query],
                                  wga_dict[reference]["INV"][query], wga_dict[reference]["SEQ"][query]]
                    var = ["GAP", "DUP", "BRK", "JMP", "INV", "SEQ"]
                    dominant_var = var[var_counts.index(max(var_counts))]
                    dominants.append(dominant_var)
            dominant_output.write("{ref}\t{dominants}\n".format(ref=reference, dominants="\t".join(dominants)))


if __name__ == "__main__":
    wga_dict, species_list = {}, []
    print("***** Species file parsing *****")
    species_parsing(args.species, species_list, wga_dict)
    print("***** Files with variation counts parsing *****")
    var_parsing(args.GAP, "GAP", wga_dict, species_list)
    var_parsing(args.DUP, "DUP", wga_dict, species_list)
    var_parsing(args.BRK, "BRK", wga_dict, species_list)
    var_parsing(args.JMP, "JMP", wga_dict, species_list)
    var_parsing(args.INV, "INV", wga_dict, species_list)
    var_parsing(args.SEQ, "SEQ", wga_dict, species_list)
    print("***** Output files creating *****")
    write_summary(wga_dict, species_list, args.out)

