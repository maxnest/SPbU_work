try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--hits', type=argparse.FileType('r'), required=True)
parser.add_argument('--hits_type', type=str, required=True,
                    help="Two type available: bbh (best BLAST hits) and leapfrog")
parser.add_argument('--output', type=str)
args = parser.parse_args()


def fasta_parsing(fasta_file, fasta_dict):
    fasta_seqs = SeqIO.parse(fasta_file, "fasta")
    for record in fasta_seqs:
        fasta_dict[record.id] = record.description.split(" ")[3].split(":")[1]


def append_query(query, hit, hits_dict):
    hit_ID = hit.split("|")[0]
    if hit_ID not in hits_dict.keys():
        hits_dict[hit_ID] = []
    hits_dict[hit_ID].append(query)


def hits_parsing(hits, hits_type, hits_dict):
    header = hits.readline()
    for line in hits:
        description = line.strip().split("\t")
        if hits_type == "leapfrog":
            query, hit = description[0], description[1]
            append_query(query, hit, hits_dict)

        elif hits_type == "bbh":
            query, hit = description[0], description[3]
            if hit != "no hits":
                append_query(query, hit, hits_dict)


def output_writing(output, fasta_dict, hits_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Query_protein_ID\tSubject_protein_ID\tSubject_gene_ID\n")
        for hit, values in hits_dict.items():
            for query in values:
                output_file.write("{query}\t{subject_prot}\t{subject_gene}\n".format(
                    query=query, subject_prot=hit, subject_gene=fasta_dict[hit]))


if __name__ == "__main__":
    fasta_dict, hits_dict = {}, {}
    fasta_parsing(args.fasta, fasta_dict)
    hits_parsing(args.hits, args.hits_type, hits_dict)
    output_writing(args.output, fasta_dict, hits_dict)