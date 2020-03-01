try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with results of sequence annotation with eggNOG-mapper")
parser.add_argument('--tab_with_selected', type=argparse.FileType('r'), required=True,
                    help="Table with sequences with normalized expression levels above threshold. "
                         "For instance, which have expression level higher than 1 scaledTPM")
parser.add_argument('--prot_2_gene_map', type=argparse.FileType('r'), required=True,
                    help="Table with 3 columns with IDs of genes, transcripts and proteins in this order")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def prot_2_gene_map_parsing(prot_2_gene_map, gene_dict):
    header = prot_2_gene_map.readline()
    for line in prot_2_gene_map:
        description = line.strip().split("\t")
        gene, transcript, protein = description[0], description[1], description[2]
        gene_dict[gene] = protein


def selected_sequences_table(table, gene_dict, selected):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene_ID, scaledTPM = description[0], description[1:]
        if gene_ID in gene_dict.keys():
            selected.append(gene_ID)
    print("{num} of sequences are selected".format(num=len(selected)))


def eggNOG_mapper_parsing(protein_dict, eggNOG):
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein_ID, annotation = description[0], description[1:]
            if len(annotation) >= 5 and len(annotation[5]) != 0:
                protein_dict[protein_ID] = annotation[5]
    print("{num} of sequences have GO-terms".format(num=len(protein_dict.keys())))


def output_creating(gene_dict, protein_dict, selected, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        for gene in selected:
            if gene_dict[gene] in protein_dict.keys():
                output.write("{gene}\t{terms}\n".format(gene=gene, terms=protein_dict[gene_dict[gene]]))


if __name__ == "__main__":
    gene_dict, protein_dict, selected = {}, {}, []
    prot_2_gene_map_parsing(args.prot_2_gene_map, gene_dict)
    selected_sequences_table(args.tab_with_selected, gene_dict, selected)
    eggNOG_mapper_parsing(protein_dict, args.eggnog)
    output_creating(gene_dict, protein_dict, selected, args.out)
