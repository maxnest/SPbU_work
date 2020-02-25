try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--gene_trans_map', type=argparse.FileType('r'), required=True,
                    help="File with 'gene_trans_map' description created by Trinity")
parser.add_argument('--clustered_fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta file with results of clusterization performed with CDHIT-est software")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def fasta_parsing(list_of_clustered, fasta):
    contigs_fasta = SeqIO.parse(fasta, 'fasta')
    for seq in contigs_fasta:
        list_of_clustered.append(seq.id.strip())


def map_parsing(map_dict, gene_trans_map):
    for line in gene_trans_map:
        description = line.strip().split('\t')
        gene, trans = description[0], description[1]
        if gene not in map_dict.keys():
            map_dict[gene] = []
            map_dict[gene].append(trans)
        else:
            map_dict[gene].append(trans)


def output_writing(map_dict, list_of_clusted, out):
    with open("{out}.trans_gene_map".format(out=out), 'a') as output:
        output.write("Transcripts\tGenes\n")
        for gene, trans_list in map_dict.items():
            for trans in trans_list:
                if trans in list_of_clusted:
                    output.write("{trans}\t{gene}\n".format(trans=trans, gene=gene))


if __name__ == "__main__":
    list_of_clustered, map_dict = [], {}
    fasta_parsing(list_of_clustered, args.clustered_fasta)
    map_parsing(map_dict, args.gene_trans_map)
    output_writing(map_dict, list_of_clustered, args.out)

