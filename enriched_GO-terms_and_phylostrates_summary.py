try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--goea_results_merged', type=argparse.FileType('r'), required=True,
                    help="Table with merged results of GeneOntology enrichment analysis (GOEA) for phylostrata."
                         "The last column should contain information about species")
parser.add_argument('--species', type=argparse.FileType('r'), required=True,
                    help="Text file with species tags (one per line) included in table analyzed")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, species_list, table_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        goid, term, species = description[0], description[1], description[-1][1:-1]
        go_key = "{id}|{term}".format(id=goid, term=term)
        if go_key not in table_dict.keys():
            table_dict[go_key] = {species: [] for species in species_list}
        table_dict[go_key][species].append("*")


def output_writing(output, species_list, table_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("GO_ID\tGO_term\t{species}\n".format(species="\t".join(species_list)))
        for go_key, values in table_dict.items():
            species_values = [values[species][0] for species in species_list]
            output_file.write("{id}\t{term}\t{species}\n".format(
                id=go_key.split("|")[0], term=go_key.split("|")[1],
                species="\t".join(species_values)))


if __name__ == "__main__":
    species_list, table_dict = [], {}
    for line in args.species:
        species_list.append(line.strip())

    table_parsing(args.goea_results_merged, species_list, table_dict)

    for go_key, values in table_dict.items():
        for species in species_list:
            if len(values[species]) == 0:
                values[species].append("-")

    output_writing(args.output, species_list, table_dict)




