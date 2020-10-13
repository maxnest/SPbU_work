try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import collections

parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="Table with orthogroups. All empty spaces should be replaced with '-'")
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


def orthogroups_parsing(ortho, ortho_dict):
    header = ortho.readline()
    for line in ortho:
        description = line.strip().split("\t")
        orthogroup, protein_sets = description[0], description[1:]

        fgig, pvit, treg, tszi, csin, fhep, psim, \
        mlig, ofel, oviv, shae, sjap, sman, smed = protein_sets[0], protein_sets[1], protein_sets[2], \
                                                   protein_sets[3], protein_sets[4], protein_sets[5], \
                                                   protein_sets[6], protein_sets[7], protein_sets[8], \
                                                   protein_sets[9], protein_sets[10], protein_sets[11], \
                                                   protein_sets[12], protein_sets[13]

        ortho_dict[orthogroup] = {"Csinensis": csin.split(", "), "Fgigantica": fgig.split(", "),
                                  "Fhepatica": fhep.split(", "), "Mlignano": mlig.split(", "),
                                  "Ofelineus": ofel.split(", "), "Oviverrini": oviv.split(", "),
                                  "Psimillimum": psim.split(", "), "Shaematobium": shae.split(", "),
                                  "Sjaponicum": sjap.split(", "), "Smansoni": sman.split(", "),
                                  "Smediterranea": smed.split(", "), "Tregenti": treg.split(", "),
                                  "Tszidati": tszi.split(", "), "Pvittatus": pvit.split(", ")}


def levels_parsing(phylostrata, levels_dict):
    for line in phylostrata:
        description = line.strip().split("\t")
        levels_dict[description[0][1:-1]] = description[1]


def levels_in_orthogroups(ortho_dict, phylostrata_dict, levels_dict, ortho_with_levels_dict):
    levels = [level for level in levels_dict.values()]
    for orthogroup, species_values in ortho_dict.items():
        ortho_with_levels_dict[orthogroup] = {level: [] for level in set(levels)}
        for species, proteins in species_values.items():
            for protein in proteins:
                if protein != "-":
                    level = levels_dict[phylostrata_dict[species][protein]]
                    # print(level)
                    ortho_with_levels_dict[orthogroup][level].append(species)

    for orthogroup, levels_values in ortho_with_levels_dict.items():
        for level, values in levels_values.items():
            if len(values) == 0:
                values.append("-")


def patterns_comparison(ortho_dict, ortho_with_levels_dict, patterns_comparison_dict):
    species = ["Csinensis", "Fgigantica", "Fhepatica", "Mlignano", "Ofelineus", "Oviverrini",
               "Psimillimum", "Shaematobium", "Sjaponicum", "Smansoni", "Smediterranea",
               "Tregenti", "Tszidati", "Pvittatus"]

    all_common_orthogroups_for_pairs = {}
    for first_species in species:
        for second_species in species:
            pair_key = "{first}_vs_{second}".format(first=first_species, second=second_species)
            all_common_orthogroups_for_pairs[pair_key] = []
            for orthogroup, species_values in ortho_dict.items():
                if species_values[first_species][0] != "-" and species_values[second_species][0] != "-":
                    all_common_orthogroups_for_pairs[pair_key].append(orthogroup)

    for pair_key, common_orthogroups in all_common_orthogroups_for_pairs.items():
        first_species, second_species = pair_key.split("_vs_")[0], pair_key.split("_vs_")[1]
        orthogroups_with_similar_pattern = []
        for orthogroup in common_orthogroups:
            first_species_levels, second_species_levels = [], []
            for level, species_list in ortho_with_levels_dict[orthogroup].items():
                if first_species in species_list:
                    first_species_levels.append(level)

                if second_species in species_list:
                    second_species_levels.append(level)
            if collections.Counter(first_species_levels) == collections.Counter(second_species_levels):
                # Lists of phylostratigraphic levels should be identical
                orthogroups_with_similar_pattern.append(orthogroup)

        patterns_comparison_dict[pair_key] = round(len(orthogroups_with_similar_pattern)/len(common_orthogroups), 2)


def output_writing(output, ortho_with_levels_dict, patterns_comparison_dict):
    ordered_levels = ["Cellular_organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa",
                      "Bilateria", "Protostomia", "Spiralia", "Lophotrochozoa", "Platyhelminthes",
                      "Class", "Order", "Family", "Genus", "Species"]

    with open("{output}.phylostrates_in_orthogroups.tsv".format(output=output), 'a') as output_file:
        output_file.write("Orthogroup_ID\t{levels}\n".format(levels="\t".join(ordered_levels)))
        for orthogroup, levels in ortho_with_levels_dict.items():
            values = []
            for level in ordered_levels:
                if levels[level][0] != "-":
                    species_in_orthogroup = ["{species}:{count}".format(
                        species=species, count=levels[level].count(species)) for species in set(levels[level])]
                    values.append("{species_and_counts}".format(species_and_counts="|".join(species_in_orthogroup)))
                else:
                    values.append("-")
            output_file.write("{orthogroup}\t{values}\n".format(orthogroup=orthogroup, values="\t".join(values)))

    with open("{output}.orthogroups_with_one_phylostratigraphic_level.tsv".format(output=output), 'a') as one_phylostrata:
        one_phylostrata.write("Orthogroup_ID\tPhylostratigraphic_level\n")
        for orthogroup, levels in ortho_with_levels_dict.items():
            values = []
            for level, species_list in levels.items():
                if species_list[0] != "-":
                    values.append(level)
            if len(values) == 1:
                one_phylostrata.write("{orthogroup}\t{values}\n".format(orthogroup=orthogroup,
                                                                        values="|".join(set(values))))

    with open("{output}.patterns_comparison.tsv".format(output=output), 'a') as patterns_comparison:
        species = ["Csinensis", "Fgigantica", "Fhepatica", "Mlignano", "Ofelineus", "Oviverrini",
                   "Psimillimum", "Shaematobium", "Sjaponicum", "Smansoni", "Smediterranea",
                   "Tregenti", "Tszidati", "Pvittatus"]
        patterns_comparison.write("Species\t{species}\n".format(species="\t".join(species)))
        for first_species in species:
            comparison_results = []
            for second_species in species:
                pair_key = "{first}_vs_{second}".format(first=first_species, second=second_species)
                comparison_results.append(str(patterns_comparison_dict[pair_key]))
            patterns_comparison.write("{first_species}\t{comparison_results}\n".format(
                first_species=first_species, comparison_results="\t".join(comparison_results)))


if __name__ == "__main__":
    phylostrata_dict, ortho_dict, levels_dict, ortho_with_levels_dict, patterns_comparison_dict = {}, {}, {}, {}, {}
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
    orthogroups_parsing(args.ortho, ortho_dict)
    levels_parsing(args.levels, levels_dict)
    print("***** Data analysis *****")
    levels_in_orthogroups(ortho_dict, phylostrata_dict, levels_dict, ortho_with_levels_dict)
    patterns_comparison(ortho_dict, ortho_with_levels_dict, patterns_comparison_dict)
    print("***** Output writing *****")
    output_writing(args.output, ortho_with_levels_dict, patterns_comparison_dict)
