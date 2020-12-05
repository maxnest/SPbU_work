try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta-file with aminoacid sequences alignment")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()

if __name__ == "__main__":
    records = SeqIO.parse(args.fasta, "fasta")
    #   for fasta in records:
    #        print("{id}: {length}".format(id=fasta.id, length=len(fasta.seq)))
    count = SeqIO.write(records, "{output}.phylip".format(output=args.output), "phylip")
    print("Converted {count} records".format(count=count))