try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="Preticulata_summary.py output table")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(tab, class_dict, nonclass_dict):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        ID, ortho, nt, nr, swiss, domains, bio, mol, kegg, exsec, specificity = description[0], description[1], \
                                                                                description[2], description[4], \
                                                                                description[6], description[8], \
                                                                                description[9], description[10], \
                                                                                description[11], description[12], \
                                                                                description[13]
        if exsec == "Classical":
            class_dict["Contigs"].append(ID)
            if ortho != '-':
                class_dict["Ortho"].append(ID)
            if nt != "No hit":
                class_dict["Anno"].append(ID)
            protein_annotation = [nr, swiss, bio, mol, kegg]
            for hit in protein_annotation:
                if hit != "No hit" and hit != "-":
                    class_dict["Anno"].append(ID)
            if domains != '-':
                class_dict["Domains"].append(ID)
            if specificity == "Externa-specific":
                class_dict["Externa"].append(ID)
            elif specificity == "Growing_stolon-specific":
                class_dict["Growing"].append(ID)
            elif specificity == "Middle_stolon_part-specific":
                class_dict["Middle"].append(ID)
            elif specificity == "Terminal_stolon_part-specific":
                class_dict["Terminal"].append(ID)
            elif specificity == "Whole_body-specific":
                class_dict["Whole"].append(ID)

        elif exsec == "Non-classical":
            nonclass_dict["Contigs"].append(ID)
            if ortho != '-':
                nonclass_dict["Ortho"].append(ID)
            if nt != "No hit":
                nonclass_dict["Anno"].append(ID)
            protein_annotation = [nr, swiss, bio, mol, kegg]
            for hit in protein_annotation:
                if hit != "No hit" and hit != "-":
                    nonclass_dict["Anno"].append(ID)
            if domains != '-':
                nonclass_dict["Domains"].append(ID)
            if specificity == "Externa-specific":
                nonclass_dict["Externa"].append(ID)
            elif specificity == "Growing_stolon-specific":
                nonclass_dict["Growing"].append(ID)
            elif specificity == "Middle_stolon_part-specific":
                nonclass_dict["Middle"].append(ID)
            elif specificity == "Terminal_stolon_part-specific":
                nonclass_dict["Terminal"].append(ID)
            elif specificity == "Whole_body-specific":
                nonclass_dict["Whole"].append(ID)


def write_output(dict, tag, out):
    with open("{out}.{tag}.tab".format(out=out, tag=tag), 'a') as output:
        output.write("Count: {contigs}\nAnno: {anno}\nOrtho: {ortho}\nDomains: {domains}\nWhole: {whole}\n"
                     "Externa: {externa}\nTerminal: {terminal}\nMiddle: {middle}\nGrowing: {growing}\n".format(
                      contigs=len(dict["Contigs"]), anno=len(set(dict["Anno"])), ortho=len(dict["Ortho"]),
                      domains=len(dict["Domains"]), whole=len(dict["Whole"]), externa=len(dict["Externa"]),
                      terminal=len(dict["Terminal"]), middle=len(dict["Middle"]), growing=len(dict["Growing"])))


if __name__ == "__main__":
    class_dict, nonclass_dict = {"Contigs": [], "Anno": [], "Ortho": [], "Domains": [],
                                 "Whole": [], "Externa": [], "Terminal": [], "Middle": [], "Growing": []}, \
                                {"Contigs": [], "Anno": [], "Ortho": [], "Domains": [],
                                 "Whole": [], "Externa": [], "Terminal": [], "Middle": [], "Growing": []}
    print("***** Input table parsing *****")
    table_parsing(args.tab, class_dict, nonclass_dict)
    print("***** Output files creating *****")
    write_output(class_dict, "potential_classical_exsec", args.out)
    write_output(nonclass_dict, "potential_nonclassical_exsec", args.out)