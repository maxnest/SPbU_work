try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str)
args = parser.parse_args()


def fasta_dict(file, dict):
    fasta_seqs = SeqIO.parse(file, "fasta")
    for fasta in fasta_seqs:
        name, sequence = fasta.id.split(" ")[0], fasta.seq
        dict[name] = {"seq": sequence, "hit_name": [], "hit_description": []}


def append_annotation(dict, tab):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        contig_ID, hit_name = description[0], description[3]
        dict[contig_ID]["hit_name"].extend(hit_name.split(" "))
        if hit_name != "no hits":
            hit_description = description[4].split(" ")
            dict[contig_ID]["hit_description"].extend(hit_description)

    for id, values in dict.items():
        if len(values["hit_description"]) == 0:
            values["hit_description"].extend(["no", "description"])


def write_output(dict, tag):
    with open("{output}.fasta".format(output=tag), 'a') as output:
        for id, values in dict.items():
            output.write(">{id}_{name}_{description}\n{seq}\n".format(id=id, name="_".join(values["hit_name"]),
                                                                      description="_".join(values["hit_description"]),
                                                                      seq=values["seq"]))


if __name__ == "__main__":
    contig_dict = {}
    print("***** Parsing fasta *****")
    fasta_dict(args.fasta, contig_dict)
    print("***** Append annotation to contigs *****")
    append_annotation(contig_dict, args.tab)
    print("***** Output file creating *****")
    write_output(contig_dict, args.output)
