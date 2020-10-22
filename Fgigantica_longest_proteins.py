try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--proteins', type=argparse.FileType('r'), required=True,
                    help="Fasta file with all aminoacid sequences for Fasciola gigantica")
args = parser.parse_args()


def fasta_parsing(fasta, fasta_dict):
    fasta_file = SeqIO.parse(fasta, "fasta")
    for el in fasta_file:
        seq_name, sequence = el.id, el.seq
        cluster_ID = seq_name.split(";")[0]
        if cluster_ID not in fasta_dict.keys():
            fasta_dict[cluster_ID] = {}
        fasta_dict[cluster_ID][seq_name] = {"sequence": sequence, "length": len(sequence)}


def longest_proteins_searching(fasta_dict, longest_sequences):
    for cluster_ID, seq_values in fasta_dict.items():
        all_lengths = [values["length"] for values in seq_values.values()]

        for seq_name, values in seq_values.items():
            if values["length"] == max(all_lengths):
                longest_sequences.append(seq_name)
                break


def output_writing(fasta_dict, longest_sequences):
    with open("Fgigantica_longest_pep_per_cluster.fasta", 'a') as output_file:
        for cluster_ID, seq_values in fasta_dict.items():
            for seq_name, values in seq_values.items():
                if seq_name in longest_sequences:
                    output_file.write(">{id}\n{sequences}\n".format(id=seq_name, sequences=values["sequence"]))


if __name__ == "__main__":
    fasta_dict, longest_sequences = {}, []
    fasta_parsing(args.proteins, fasta_dict)
    longest_proteins_searching(fasta_dict, longest_sequences)
    output_writing(fasta_dict, longest_sequences)