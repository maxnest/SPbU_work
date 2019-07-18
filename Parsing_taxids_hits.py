try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--tag', type=str, required=True)
args = parser.parse_args()

hits_dict = {}


def parsing(table, dict):
    for line in table:
        description = line.strip().split("\t")
        prot_ID, hit_desc = description[0], description[2].split("|")
        if len(hit_desc) != 0:
            if "Eukaryota" in hit_desc:
                dict[prot_ID] = hit_desc[4]
            else:
                dict[prot_ID] = hit_desc[2]
        else:
            dict[prot_ID] = "No_hit"
            print(hits_dict)

        # taxids_export.txt:
        # root|cellular organisms|Eukaryota|Opisthokonta|Fungi|
        # root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|
        # root|cellular organisms|Bacteria|FCB group|Bacteroidetes/Chlorobi group|


def write_output(dict, tag):
    with open("{tag}.tsv".format(tag=tag), 'a') as table:
        table.write("Protein_ID\tHit\n")
        for key, value in dict.items():
            table.write("{prot}\t{hit}\n".format(prot=key, hit=value))


if __name__ == "__main__":
    parsing(args.tab, hits_dict)
    write_output(hits_dict, args.tag)