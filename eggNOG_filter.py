try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
                    help="Table with results of sequence annotation performed with help of eggNOG-mapper")
parser.add_argument('--exclude', type=argparse.FileType('r'), required=True,
                    help="Text file with list of taxons (one per line) that should be excluded")
parser.add_argument('--map', type=argparse.FileType('r'), required=True,
                    help="Table with map of IDs: gene transcript protein")
parser.add_argument('--exp', type=argparse.FileType('r'), required=True,
                    help="Table with expression values")
parser.add_argument('--transcripts', type=argparse.FileType('r'), required=True,
                    help="Fasta file with nucleotide sequences")
parser.add_argument('--proteins', type=argparse.FileType('r'), required=True,
                    help="Fasta file with aminoacid sequences")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def eggNOG_mapper_parsing(eggnog, eggnog_dict):
    for line in eggnog:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein, taxon = description[0], description[17]
            eggnog_dict[protein] = taxon


def exclude_parsing(exclude, list_of_excluded):
    for line in exclude:
        description = line.strip().split("\t")
        list_of_excluded.append(description[0])


def map_parsing_and_filtering(map, list_of_excluded, eggnog_dict, map_dict, output):
    header = map.readline().strip().split("\t")
    for line in map:
        description = line.strip().split("\t")
        gene, transcript, protein = description[0], description[1], description[2]
        if protein not in eggnog_dict.keys() or eggnog_dict[protein] not in list_of_excluded:
            map_dict[gene] = {"transcript": transcript, "protein": protein}

    with open("{output}.eggnog_filter.gene_map.tsv".format(output=output), 'a') as output_map:
        output_map.write("{header}\n".format(header="\t".join(header)))
        for gene, products in map_dict.items():
            output_map.write("{gene}\t{transcript}\t{protein}\n".format(gene=gene, transcript=products["transcript"],
                                                                        protein=products["protein"]))


def expression_parsing_and_filtering(exp, map_dict, output):
    header = exp.readline().strip().split("\t")
    with open("{output}.eggnog_filter.averaged_scaledTPM.tsv".format(output=output), 'a') as output_exp:
        output_exp.write("{header}\n".format(header="\t".join(header)))
        for line in exp:
            description = line.strip().split("\t")
            gene, expression = description[0], description[1:]
            if gene in map_dict.keys():
                output_exp.write("{gene}\t{expression}\n".format(gene=gene, expression="\t".join(expression)))


def fasta_parser(fasta, type, map_dict, output):
    transcripts = [map_dict[gene]["transcript"] for gene in map_dict.keys()]
    proteins = [map_dict[gene]["protein"] for gene in map_dict.keys()]
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    with open("{output}.eggnog_filter.{type}.fasta".format(output=output, type=type), 'a') as output_fasta:
        for fasta_seq in fasta_seqs:
            name, sequence = fasta_seq.id, fasta_seq.seq
            if name in transcripts or name in proteins:
                output_fasta.write(">{name}\n{seq}\n".format(name=name, seq=sequence))


if __name__ == "__main__":
    eggnog_dict, list_of_excluded, map_dict = {}, [], {}
    print("***** eggNOG table parsing *****")
    eggNOG_mapper_parsing(args.eggnog, eggnog_dict)
    exclude_parsing(args.exclude, list_of_excluded)
    print("These taxons will be excluded: {excluded}".format(excluded="\t".join(list_of_excluded)))
    print("***** Filter applying *****")
    map_parsing_and_filtering(args.map, list_of_excluded, eggnog_dict, map_dict, args.output)
    expression_parsing_and_filtering(args.exp, map_dict, args.output)
    fasta_parser(args.transcripts, "nucl", map_dict, args.output)
    fasta_parser(args.proteins, "prot", map_dict, args.output)
