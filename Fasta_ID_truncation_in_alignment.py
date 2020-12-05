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
    fasta_dict = {}
    records = SeqIO.parse(args.fasta, "fasta")
    for fasta in records:
        fasta_dict[fasta.id] = {"fasta_new_id": "{genus}{species}".format(genus=fasta.id.split("_")[0][0],
                                                                          species=fasta.id.split("_")[1][0:9]),
                                "fasta_seq": fasta.seq}
    with open("{output}.fasta".format(output=args.output), 'a') as output_fasta:
        print("Renaming: ")
        for fasta_id, values in fasta_dict.items():
            print("{old} --> {new}".format(old=fasta_id, new=values["fasta_new_id"]))
            output_fasta.write(">{new}\n{seq}\n".format(new=values["fasta_new_id"], seq=values["fasta_seq"]))

