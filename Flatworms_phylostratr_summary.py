try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
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
parser.add_argument('--levels', type=argparse.FileType('r'), required=True,
                    help="Table with description of phylostratigraphic levels")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def phylostratr_table_parsing(table, species_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        if mrca_name not in species_dict.keys():
            species_dict[mrca_name] = []
        species_dict[mrca_name].append(protein_ID)


def levels_parsing(levels, levels_dict):
    for line in levels:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def merge_results(merged_dict, levels, ordered_levels, csinensis_dict, fgigantica_dict, fhepatica_dict,
                  mlignano_dict, ofelineus_dict, oviverrini_dict, psimillimum_dict, pvittatus_dict,
                  shaematobium_dict, sjaponicum_dict, smansoni_dict, smediterranea_dict,
                  tregenti_dict, tszidati_dict):
    species_dicts = {"Csinensis": csinensis_dict, "Fgigantica": fgigantica_dict, "Fhepatica": fhepatica_dict,
                     "Mlignano": mlignano_dict, "Ofelineus": ofelineus_dict, "Oviverrini": oviverrini_dict,
                     "Psimillimum": psimillimum_dict, "Pvittatus": pvittatus_dict,
                     "Shaematobium": shaematobium_dict, "Sjaponicum": sjaponicum_dict, "Smansoni": smansoni_dict,
                     "Smediterranea": smediterranea_dict, "Tregenti": tregenti_dict, "Tszidati": tszidati_dict}

    for level in ordered_levels:
        merged_dict[level] = {"Csinensis": [], "Fgigantica": [], "Fhepatica": [], "Mlignano": [],
                              "Ofelineus": [], "Oviverrini": [], "Psimillimum": [], "Pvittatus": [],
                              "Shaematobium": [], "Sjaponicum": [], "Smansoni": [], "Smediterranea": [],
                              "Tregenti": [], "Tszidati": []}

    for species, species_dict in species_dicts.items():
        protein_number = 0
        for mrca_name, proteins in species_dict.items():
            merged_dict[levels[mrca_name]][species].append(len(proteins))
            protein_number += len(proteins)
        # append per cents
        for level in merged_dict.keys():
            if len(merged_dict[level][species]) != 0:
                merged_dict[level][species].append((int(merged_dict[level][species][0])/protein_number)*100)
            else:
                merged_dict[level][species].extend([0, 0])


def write_output(output, ordered_levels, merged_dict):
    with open("{output}.phylostratr_summary.tsv".format(output=output), 'a') as phylostratr_summary:
        phylostratr_summary.write("Levels\tCsinensis\tFgigantica\tFhepatica\tMlignano\tOfelineus\t"
                                  "Oviverrini\tPsimillimum\tPvittatus\tShaematobium\tSjaponicum\tSmansoni\t"
                                  "Smediterranea\tTregenti\tTszidati\n")
        for level in ordered_levels:
            values = ["{num}({percent}%)".format(num=merged_dict[level][species][0],
                                                 percent=round(merged_dict[level][species][1], 2)) for species in
                      ["Csinensis", "Fgigantica", "Fhepatica", "Mlignano", "Ofelineus", "Oviverrini", "Psimillimum",
                       "Pvittatus", "Shaematobium", "Sjaponicum", "Smansoni", "Smediterranea", "Tregenti", "Tszidati"]]
            phylostratr_summary.write("{level}\t{values}\n".format(level=level, values="\t".join(values)))


if __name__ == "__main__":
    ordered_levels = ["Cellular_organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria",
                      "Protostomia", "Spiralia", "Lophotrochozoa", "Platyhelminthes", "Class",
                      "Order", "Family", "Genus", "Species"]
    csinensis_dict, fgigantica_dict, fhepatica_dict, mlignano_dict, ofelineus_dict, oviverrini_dict, psimillimum_dict, \
    pvittatus_dict, shaematobium_dict, sjaponicum_dict, smansoni_dict, smediterranea_dict, tregenti_dict, \
    tszidati_dict = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    levels_dict, merged_dict = {}, {}
    print("***** Input files parsing *****")
    phylostratr_table_parsing(args.csinensis, csinensis_dict)
    phylostratr_table_parsing(args.fgigantica, fgigantica_dict)
    phylostratr_table_parsing(args.fhepatica, fhepatica_dict)
    phylostratr_table_parsing(args.mlignano, mlignano_dict)
    phylostratr_table_parsing(args.ofelineus, ofelineus_dict)
    phylostratr_table_parsing(args.oviverrini, oviverrini_dict)
    phylostratr_table_parsing(args.psimillimum, psimillimum_dict)
    phylostratr_table_parsing(args.pvittatus, pvittatus_dict)
    phylostratr_table_parsing(args.shaematobium, shaematobium_dict)
    phylostratr_table_parsing(args.sjaponicum, sjaponicum_dict)
    phylostratr_table_parsing(args.smansoni, smansoni_dict)
    phylostratr_table_parsing(args.smediterranea, smediterranea_dict)
    phylostratr_table_parsing(args.tregenti, tregenti_dict)
    phylostratr_table_parsing(args.tszidati, tszidati_dict)
    levels_parsing(args.levels, levels_dict)
    print("***** Results merging *****")
    merge_results(merged_dict, levels_dict, ordered_levels, csinensis_dict, fgigantica_dict, fhepatica_dict,
                  mlignano_dict, ofelineus_dict, oviverrini_dict, psimillimum_dict, pvittatus_dict,
                  shaematobium_dict, sjaponicum_dict, smansoni_dict, smediterranea_dict, tregenti_dict, tszidati_dict)
    print("***** Output creating *****")
    write_output(args.output, ordered_levels, merged_dict)
