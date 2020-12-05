try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--domains', type=argparse.FileType('r'), required=True,
                    help="Text file with wanted domains (one per line)")
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta-file with aminoacid sequences")
parser.add_argument('--species_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def text_file_parsing(domains, list_of_domains):
    for line in domains:
        list_of_domains.append(line.strip())


def sequences_ordering(fasta, list_of_domains, ordered_seqs):
    fasta_seqs_dict = {}

    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for fasta in fasta_seqs:
        fasta_seqs_dict[fasta.id.split("_")[0]] = str(fasta.seq)

    for domain_ID in list_of_domains:
        if domain_ID in fasta_seqs_dict.keys():
            ordered_seqs.append(fasta_seqs_dict[domain_ID])


def output_writing(output, species_tag, ordered_seqs):
    with open("{output}.fasta".format(output=output), 'a') as output_fasta:
        output_fasta.write(">{species}\n{seq}\n".format(species=species_tag, seq="".join(ordered_seqs)))

    print("The total length of the final sequence: {length}".format(length=len("".join(ordered_seqs))))


if __name__ == "__main__":
    list_of_domains, ordered_seqs = [], []
    text_file_parsing(args.domains, list_of_domains)
    sequences_ordering(args.fasta, list_of_domains, ordered_seqs)
    output_writing(args.output, args.species_tag, ordered_seqs)