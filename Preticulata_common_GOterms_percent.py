try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--externa', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for Externa")
parser.add_argument('--growing', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for Growing stolon part")
parser.add_argument('--middle', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for Middle stolon part")
parser.add_argument('--terminal', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for Terminal stolon part")
parser.add_argument('--whole_body', type=argparse.FileType('r'), required=True,
                    help="Table with enriched GOterms for Whole body")
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


def unite_dict(united_dict, externa_dict, growing_dict, middle_dict, terminal_dict, whole_body_dict):
    all_terms, site_dicts = [], [externa_dict, growing_dict, middle_dict, terminal_dict, whole_body_dict]
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
                "Externa": number_of_significant[0], "Growing": number_of_significant[1],
                "Middle": number_of_significant[2], "Terminal": number_of_significant[3],
                "Whole_body": number_of_significant[4]
            }


def output_writing(output, united_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("GOterm\tExterna\tGrowing\tMiddle\tTerminal\tWhole_body\n")
        for key, values in united_dict.items():
            output_file.write("{term}\t{externa}\t{growing}\t{middle}\t{terminal}\t{whole}\n".format(
                term=key, externa=values["Externa"], growing=values["Growing"], middle=values["Middle"],
                terminal=values["Terminal"], whole=values["Whole_body"]
            ))


if __name__ == "__main__":
    externa_dict, growing_dict, middle_dict, terminal_dict, whole_body_dict = {}, {}, {}, {}, {}
    united_dict = {}
    table_parsing(args.externa, externa_dict)
    table_parsing(args.growing, growing_dict)
    table_parsing(args.middle, middle_dict)
    table_parsing(args.terminal, terminal_dict)
    table_parsing(args.whole_body, whole_body_dict)
    unite_dict(united_dict, externa_dict, growing_dict, middle_dict, terminal_dict, whole_body_dict)
    output_writing(args.output, united_dict)