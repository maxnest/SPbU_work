try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(tab, COG_dict):
    COG_keys = ["D:Cell cycle control, cell division, chromosome partitioning",
                "M:Cell wall/membrane/envelope biogenesis",
                "N:Cell motility",
                "O:Post-translational modification, protein turnover, and chaperones",
                "T:Signal transduction mechanisms",
                "U:Intracellular trafficking, secretion, and vesicular transport",
                "V:Defense mechanisms",
                "W:Extracellular structures",
                "Y:Nuclear structure",
                "Z:Cytoskeleton",
                "A:RNA processing and modification",
                "B:Chromatin structure and dynamics",
                "J:Translation, ribosomal structure and biogenesis",
                "K:Transcription",
                "L:Replication, recombination and repair",
                "C:Energy production and conversion",
                "E:Amino acid transport and metabolism",
                "F:Nucleotide transport and metabolism",
                "G:Carbohydrate transport and metabolism",
                "H:Coenzyme transport and metabolism",
                "I:Lipid transport and metabolism",
                "P:Inorganic ion transport and metabolism",
                "Q:Secondary metabolites biosynthesis, transport, and catabolism",
                "R:General function prediction only",
                "S:Function unknown"]

    for COG in COG_keys:
        COG_dict[COG.split(":")[0]] = {"description": [COG.split(":")[1]],
                                       "Externa-specific": [], "Growing_stolon-specific": [],
                                       "Middle_stolon_part-specific": [], "Terminal_stolon_part-specific": [],
                                       "Whole_body-specific": []}

    for line in tab:
        if not line.startswith("#") and not line.startswith("Contig_ID"):
            description = line.strip().split("\t")
            contig_ID, COG, jongeneel = description[0], list(description[21]), description[24]
            if jongeneel != "-" and COG[0] != "-":
                for el in COG:
                    COG_dict[el][jongeneel].append(contig_ID)


def write_output(output, COG_dict):
    with open("{out}.tsv".format(out=output), 'a') as output:
        output.write("COG\tDescription\tJongeneel:Externa\tJongeneel:Growing_stolon\tJongeneel:Middle_stolon\t"
                     "Jongeneel:Terminal_stolon\tJongeneel:Whole_body\n")
        for COG, values in COG_dict.items():
            output.write("{cog}\t{description}\t{externa}\t{growing}\t{middle}\t{terminal}\t{whole}\n".format(
                cog=COG, description=values["description"][0], externa=len(values["Externa-specific"]),
                growing=len(values["Growing_stolon-specific"]), middle=len(values["Middle_stolon_part-specific"]),
                terminal=len(values["Terminal_stolon_part-specific"]), whole=len(values["Whole_body-specific"])
            ))


if __name__ == "__main__":
    COG_dict = {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, COG_dict)
    print("***** Output writing *****")
    write_output(args.output, COG_dict)