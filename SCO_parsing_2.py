try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--sco_csv', type=argparse.FileType('r'), required=True)
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta with all amino acids sequences used in OrthoFinder analysis")
args = parser.parse_args()


def orthogroups_parsing(sco_csv, orthodict):
    for line in sco_csv:
        description = line.strip().split("\t")
        group_ID, proteins_ID = description[0], description[1:]
        orthodict[group_ID] = proteins_ID


def write_fasta(orthodict, fasta_seqs):
    proteins_ID, fasta_dict = [], {}
    for key, value in orthodict.items():
        proteins_ID.extend(value)

    for fasta in fasta_seqs:
        name, sequence = fasta.id, fasta.seq
        if name in proteins_ID:
            fasta_dict[name] = sequence

    for key, value in ortodict.items():
        with open("{key}.fasta".format(key=key), 'a') as output:
            output.write(">Csin\n{seq}\n".format(seq=fasta_dict[value[0]]))
            output.write(">Fhep\n{seq}\n".format(seq=fasta_dict[value[1]]))
            output.write(">Mlig\n{seq}\n".format(seq=fasta_dict[value[2]]))
            output.write(">Oviv\n{seq}\n".format(seq=fasta_dict[value[3]]))
            output.write(">Shae\n{seq}\n".format(seq=fasta_dict[value[4]]))
            output.write(">Sjap\n{seq}\n".format(seq=fasta_dict[value[5]]))
            output.write(">Sman\n{seq}\n".format(seq=fasta_dict[value[6]]))
            output.write(">Smed\n{seq}\n".format(seq=fasta_dict[value[7]]))


if __name__ == "__main__":
    fasta_seqs = SeqIO.parse(args.fasta, "fasta")
    ortodict = {}
    print("*** Parsing SinlgeCopyOrthogroups ***")
    orthogroups_parsing(args.sco_csv, ortodict)
    print("*** Output fasta creating ***")
    write_fasta(ortodict, fasta_seqs)
