try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--ref', type=argparse.FileType('r'), required=True,
                    help="Fasta file with reference sequences")
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="Text file with structural variation description: dnadiff AND fgrep '***SV***'")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def fasta_reading(fasta, fasta_dict):
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for fasta in fasta_seqs:
        name, sequence = fasta.id, fasta.seq
        fasta_dict[name] = sequence


def sv_extraction(tab, fasta_dict, sv_dict):
    for line in tab:
        description = line.strip().split("\t")
        ID, start, end, length = description[0], int(description[2]), int(description[3]), int(description[4])
        if len(fasta_dict[ID][start - 1:end]) == length:
            sv_dict["{ID}_fr_{start}_to_{end}".format(ID=ID, start=start, end=end)] = fasta_dict[ID][start - 1:end]


def write_output(sv_dict, out):
    with open("{out}.fasta".format(out=out), 'a') as output:
        for key, value in sv_dict.items():
            output.write(">{name}\n{seq}\n".format(name=key, seq=value))


if __name__ == "__main__":
    fasta_dict, sv_dict = {}, {}
    print("***** Fasta file reading *****")
    fasta_reading(args.ref, fasta_dict)
    print("***** Structural variation extraction *****")
    sv_extraction(args.tab, fasta_dict, sv_dict)
    print("***** Output file creating *****")
    write_output(sv_dict, args.out)