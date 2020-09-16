try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="Table with orthogroups")
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="Table with description of phylostratigraphic levels")
parser.add_argument('--csinensis', type=argparse.FileType('r'), required=True,
                    help="Table with results of C.sinensis phylostratigraphy analysis")
parser.add_argument('--psimillimum', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum phylostratigraphy analysis")
parser.add_argument('--fhepatica', type=argparse.FileType('r'), required=True,
                    help="Table with results of F.hepatica phylostratigraphy analysis")
parser.add_argument('--fgigantica', type=argparse.FileType('r'), required=True,
                    help="Table with results of F.gigantica phylostratigraphy analysis")
parser.add_argument('--ofelineus', type=argparse.FileType('r'), required=True,
                    help="Table with results of O.felineus phylostratigraphy analysis")
parser.add_argument('--oviverrini', type=argparse.FileType('r'), required=True,
                    help="Table with results of O.viverrini phylostratigraphy analysis")
parser.add_argument('--shaematobium', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.haematobium phylostratigraphy analysis")
parser.add_argument('--sjaponicum', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.japonicum phylostratigraphy analysis")
parser.add_argument('--smansoni', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.mansoni phylostratigraphy analysis")
parser.add_argument('--tregenti', type=argparse.FileType('r'), required=True,
                    help="Table with results of T.regenti phylostratigraphy analysis")
parser.add_argument('--tszidati', type=argparse.FileType('r'), required=True,
                    help="Table with results of T.szidati phylostratigraphy analysis")
parser.add_argument('--smediterranea', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.mediterranea phylostratigraphy analysis")
parser.add_argument('--mlignano', type=argparse.FileType('r'), required=True,
                    help="Table with results of M.lignano phylostratigraphy analysis")
parser.add_argument('--pvittatus', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.vittatus phylostratigraphy analysis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_table_parsing(table, species, phylostrata_dict):
    header = table.readline()
    if species not in phylostrata_dict.keys():
        phylostrata_dict[species] = {}
    for line in table:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0], description[1], description[2], description[3]
        phylostrata_dict[species][protein_ID[1:-1]] = mrca_name[1:-1]


def orthogroups_parsing(ortho, phylostrata_dict, ortho_dict):
    header = ortho.readline()
    for line in ortho:
        description = line.strip().split("\t")
        orthogroup, protein_sets = description[0], description[1:]
        ortho_dict[orthogroup] = {"Csinensis": [], "Fgigantica": [], "Fhepatica": [], "Mlignano": [], "Ofelineus": [],
                                  "Oviverrini": [], "Psimillimum": [], "Shaematobium": [], "Sjaponicum": [],
                                  "Smansoni": [], "Smediterranea": [], "Tregenti": [], "Tszidati": [], "Pvittatus": []}
        for protein_set in protein_sets:
            for species, values in phylostrata_dict.items():
                all_keys = [key for key in values.keys()]
                if protein_set.split(", ")[0] in all_keys:
                    ortho_dict[orthogroup][species].extend(protein_set.split(", "))

        for species, proteins in ortho_dict[orthogroup].items():
            if len(proteins) == 0:
                proteins.append("-")


def levels_parsing(phylostrata, levels_dict):
    for line in phylostrata:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def phylostrates_in_orthogroups(ortho_dict, phylostrata_dict, levels_dict, ortho_with_phylostrates_dict):
    for orthogroup, values in ortho_dict.items():
        ortho_with_phylostrates_dict[orthogroup] = {"Csinensis": [], "Fgigantica": [], "Fhepatica": [],
                                                    "Mlignano": [], "Ofelineus": [], "Oviverrini": [],
                                                    "Psimillimum": [], "Shaematobium": [], "Sjaponicum": [],
                                                    "Smansoni": [], "Smediterranea": [], "Tregenti": [],
                                                    "Tszidati": [], "Pvittatus": []}
        for species, proteins in values.items():
            if proteins[0] != "-":
                for protein in proteins:
                    if protein in phylostrata_dict[species].keys():
                        ortho_with_phylostrates_dict[orthogroup][species].append(
                            levels_dict[phylostrata_dict[species][protein]])
            else:
                ortho_with_phylostrates_dict[orthogroup][species].append("-")


def output_writing(ortho_with_phylostrates_dict, output):
    with open("{output}.phylostrates_in_orthogroups.tsv".format(output=output), 'a') as ortho_output:
        ortho_output.write("Orthogroup\tCsinensis\tFgigantica\tFhepatica\tMlignano\tOfelineus\tOviverrini\t"
                           "Psimillimum\tShaematobium\tSjaponicum\tSmansoni\tSmediterranea\tTregenti\tTszidati\t"
                           "Pvittatus\n")
        for orthogroup, values in ortho_with_phylostrates_dict.items():
            ortho_output.write("{ortho}\t{csin}\t{fgig}\t{fhep}\t{mlig}\t{ofel}\t{oviv}\t{psim}\t{shae}\t{sjap}\t"
                               "{sman}\t{smed}\t{treg}\t{tszi}\t{pvit}\n".format(
                                ortho=orthogroup,
                                csin=",".join(["{level}:{count}".format(level=level, count=values["Csinensis"].count(level))
                                    for level in set(values["Csinensis"]) if level != '-']),

                                fgig=",".join(["{level}:{count}".format(level=level, count=values["Fgigantica"].count(level))
                                    for level in set(values["Fgigantica"]) if level != '-']),

                                fhep=",".join(["{level}:{count}".format(level=level, count=values["Fhepatica"].count(level))
                                    for level in set(values["Fhepatica"]) if level != '-']),

                                mlig=",".join(["{level}:{count}".format(level=level, count=values["Mlignano"].count(level))
                                    for level in set(values["Mlignano"]) if level != '-']),

                                ofel=",".join(["{level}:{count}".format(level=level, count=values["Ofelineus"].count(level))
                                    for level in set(values["Ofelineus"]) if level != '-']),

                                oviv=",".join(["{level}:{count}".format(level=level, count=values["Oviverrini"].count(level))
                                    for level in set(values["Oviverrini"]) if level != '-']),

                                psim=",".join(["{level}:{count}".format(level=level, count=values["Psimillimum"].count(level))
                                    for level in set(values["Psimillimum"]) if level != '-']),

                                shae=",".join(["{level}:{count}".format(level=level, count=values["Shaematobium"].count(level))
                                    for level in set(values["Shaematobium"]) if level != '-']),

                                sjap=",".join(["{level}:{count}".format(level=level, count=values["Sjaponicum"].count(level))
                                    for level in set(values["Sjaponicum"]) if level != '-']),

                                sman=",".join(["{level}:{count}".format(level=level, count=values["Smansoni"].count(level))
                                    for level in set(values["Smansoni"]) if level != '-']),

                                smed=",".join(["{level}:{count}".format(level=level, count=values["Smediterranea"].count(level))
                                    for level in set(values["Smediterranea"]) if level != '-']),

                                treg=",".join(["{level}:{count}".format(level=level, count=values["Tregenti"].count(level))
                                    for level in set(values["Tregenti"]) if level != '-']),

                                tszi=",".join(["{level}:{count}".format(level=level, count=values["Tszidati"].count(level))
                                    for level in set(values["Tszidati"]) if level != '-']),

                                pvit=",".join(["{level}:{count}".format(level=level, count=values["Pvittatus"].count(level))
                                    for level in set(values["Pvittatus"]) if level != '-'])))

    with open("{output}.orthogroups_with_one_phylostrata.tsv".format(output=output), 'a') as one_phylostrata:
        one_phylostrata.write("Orthogroup\tPhylostrata\n")
        for orthogroup, values in ortho_with_phylostrates_dict.items():
            all_phylostrates_in_group = []
            for species, phylostrates in values.items():
                all_phylostrates_in_group.extend(phylostrates)
            if "-" in all_phylostrates_in_group:
                all_phylostrates_in_group.remove("-")

            if len(set(all_phylostrates_in_group)) == 1:
                one_phylostrata.write("{ortho}\t{phylostrata}\n".format(ortho=orthogroup,
                                                                        phylostrata=all_phylostrates_in_group[0]))


if __name__ == "__main__":
    phylostrata_dict, ortho_dict, levels_dict, ortho_with_phylostrates_dict = {}, {}, {}, {}
    print("***** Phylostratr tables parsing *****")
    phylostratr_table_parsing(args.csinensis, "Csinensis", phylostrata_dict)
    phylostratr_table_parsing(args.psimillimum, "Psimillimum", phylostrata_dict)
    phylostratr_table_parsing(args.fhepatica, "Fhepatica", phylostrata_dict)
    phylostratr_table_parsing(args.fgigantica, "Fgigantica", phylostrata_dict)
    phylostratr_table_parsing(args.ofelineus, "Ofelineus", phylostrata_dict)
    phylostratr_table_parsing(args.oviverrini, "Oviverrini", phylostrata_dict)
    phylostratr_table_parsing(args.shaematobium, "Shaematobium", phylostrata_dict)
    phylostratr_table_parsing(args.sjaponicum, "Sjaponicum", phylostrata_dict)
    phylostratr_table_parsing(args.smansoni, "Smansoni", phylostrata_dict)
    phylostratr_table_parsing(args.tregenti, "Tregenti", phylostrata_dict)
    phylostratr_table_parsing(args.tszidati, "Tszidati", phylostrata_dict)
    phylostratr_table_parsing(args.smediterranea, "Smediterranea", phylostrata_dict)
    phylostratr_table_parsing(args.mlignano, "Mlignano", phylostrata_dict)
    phylostratr_table_parsing(args.pvittatus, "Pvittatus", phylostrata_dict)
    print("***** Orthogroups and levels parsing *****")
    orthogroups_parsing(args.ortho, phylostrata_dict, ortho_dict)
    levels_parsing(args.levels, levels_dict)
    print("***** Orthogroups analysis *****")
    phylostrates_in_orthogroups(ortho_dict, phylostrata_dict, levels_dict, ortho_with_phylostrates_dict)
    output_writing(ortho_with_phylostrates_dict, args.output)