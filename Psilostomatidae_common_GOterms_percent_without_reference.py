try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--redia', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for redia stage")
parser.add_argument('--cercaria', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for cercaria stage")
parser.add_argument('--marita', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for marita stage")
parser.add_argument('--species_tag', type=str, required=True, help="For instance: Psimillimum")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, go_dict):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        goid, term, annotated, significant = description[0], description[1], int(description[2]), int(description[3])
        go_dict[goid] = {"Term": term, "Annotated": annotated, "Significant": significant}


def append_number(list, key, site_dict):
    if key in site_dict.keys():
        list.append(round((site_dict[key]["Significant"]/site_dict[key]["Annotated"])*100, 2))    # NB: percent
    else:
        list.append(0)


def unite_dict(united_dict, redia_dict, cercaria_dict, marita_dict):
    all_terms, site_dicts = [], [redia_dict, cercaria_dict, marita_dict]
    for site_dict in site_dicts:
        all_terms.extend([key for key in site_dict.keys()])
    for key in set(all_terms):
        number_of_significant = []
        for site_dict in site_dicts:
            append_number(number_of_significant, key, site_dict)
        number_of_nonzero = 0
        for number in number_of_significant:
            if number != 0:
                number_of_nonzero += 1
        if number_of_nonzero > 1:
            #print(number_of_significant)
            united_dict["{go}|{term}".format(
                go=key,
                term=site_dicts[number_of_significant.index(
                    next(filter(lambda x: x != 0, number_of_significant)))][key]["Term"])] = {
                "Redia": number_of_significant[0], "Cercaria": number_of_significant[1],
                "Marita": number_of_significant[2]}


def output_writing(output, species_tag, united_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("GOterm\t{species}_redia\t{species}_cercaria\t{species}_marita\n".format(species=species_tag))
        for key, values in united_dict.items():
            output_file.write("{term}\t{redia}\t{cercaria}\t{marita}\n".format(
                term=key, redia=values["Redia"], cercaria=values["Cercaria"], marita=values["Marita"]))


if __name__ == "__main__":
    redia_dict, cercaria_dict, marita_dict = {}, {}, {}
    united_dict = {}
    table_parsing(args.redia, redia_dict)
    table_parsing(args.cercaria, cercaria_dict)
    table_parsing(args.marita, marita_dict)
    unite_dict(united_dict, redia_dict, cercaria_dict, marita_dict)
    output_writing(args.output, args.species_tag, united_dict)