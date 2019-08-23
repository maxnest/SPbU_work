try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--smes_gff3', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def fasta_parsing(fasta, fasta_dict):
    for line in fasta:
        if line.startswith(">"):
            fasta_dict[line.strip()[1:]] = {"seq": next(fasta).strip(), "gene_ID": []}


def gff_parsing(gff3, sp_dict):
    for line in gff3:
        if not line.startswith("#"):
            description = line.strip().split("\t")

            if description[2] == "gene":
                scf, start, end = description[0], int(description[3]), int(description[4])
                name = description[-1].split("=")[1]

                sp_dict[name] = {"scf": scf, "start": start, "end": end}


def protein2gene(fasta_dict, gff_dict):
    # fasta : >dd_Smes_g4_100:1134689-1136073.p1
    # gene ID: dd_Smes_g4_67	AUGUSTUS	gene	966008	971956	0.01	-	.	ID=SMESG000067732.1
    for protein, values in fasta_dict.items():
        protein_ID, number = protein.split(".")[0], protein.split(".")[1]
        description = protein_ID.split("_")
        scf = "{el1}_{el2}_{el3}_{el4}".format(el1=description[0], el2=description[1], el3=description[2],
                                               el4=description[3].split(":")[0])
        coordinates = description[3].split(":")[1]
        start, end = int(coordinates.split("-")[0]), int(coordinates.split("-")[1])
        for gene, gene_value in gff_dict.items():
            if gene_value["scf"] == scf:
                if gene_value["start"] <= start + 1 and gene_value["end"] >= end:
                    values["gene_ID"].append("{gene}.{num}".format(gene=gene, num=number))

    print(fasta_dict)

def write_fasta(fasta_dict, output):
    with open("{output}.fasta".format(output=output), 'a') as out_fasta:
        for protein, values in fasta_dict.items():
            out_fasta.write(">{gene}\n{seq}\n".format(gene=values["gene_ID"][0], seq=values["seq"]))


if __name__ == "__main__":
    fasta_dict, gff_dict = {}, {}
    print("*** Fasta parsing ***")
    fasta_parsing(args.fasta, fasta_dict)
    print("*** GFF3 parsing ***")
    gff_parsing(args.smes_gff3, gff_dict)
    print("*** Protein 2 Gene ***")
    protein2gene(fasta_dict, gff_dict)
    print("*** Output file creating ***")
    write_fasta(fasta_dict, args.output)
