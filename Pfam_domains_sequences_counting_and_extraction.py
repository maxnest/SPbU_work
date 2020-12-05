try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--domains', type=argparse.FileType('r'), required=True,
                    help="Text file with wanted domains (one per line)")
parser.add_argument('--arch', type=argparse.FileType('r'), required=True,
                    help="Table with domain architecture of the proteins")
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta-file with aminoacid sequences")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def text_file_parsing(domains, list_of_domains):
    for line in domains:
        list_of_domains.append(line.strip())


def domain_arch_parsing(arch, arch_dict):
    header = arch.readline()
    for line in arch:
        description = line.strip().split("\t")
        protein_ID, domain_name, domain_ID, start, end = description[0], description[1], \
            description[2].split(".")[0], int(description[3]), int(description[4])

        if domain_ID not in arch_dict.keys():
            arch_dict[domain_ID] = {}

        key = "{protein_ID}_from_{start}_to_{end}".format(protein_ID=protein_ID, start=start, end=end)
        arch_dict[domain_ID][key] = {"protein_ID": protein_ID, "start": start, "end": end, "sequence": ""}


def sequences_extraction(fasta, arch_dict):
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    fasta_seqs_dict = {}

    for fasta in fasta_seqs:
        fasta_seqs_dict[fasta.id] = str(fasta.seq)

    for domain_ID, protein_values in arch_dict.items():
        for key, protein_value in protein_values.items():
            protein_ID, start, end = protein_value["protein_ID"], protein_value["start"], protein_value["end"]
            protein_value["sequence"] += fasta_seqs_dict[protein_ID][start - 1: end]


def output_writing(output, arch_dict, list_of_domains):
    with open("{output}.wanted_domains_counts_and_protein_IDs.tsv".format(output=output), 'a') as output_file:
        output_file.write("Domain_IDs\tCount\tProtein_IDs_and_coordinates\n")
        for domain_ID in list_of_domains:
            if domain_ID in arch_dict.keys():
                output_file.write("{domain}\t{count}\t{proteins}\n".format(
                    domain=domain_ID, count=len(arch_dict[domain_ID].keys()),
                    proteins=";".join([key for key in arch_dict[domain_ID].keys()])
                ))

    for domain_ID in list_of_domains:
        if domain_ID in arch_dict.keys():
            with open("{output}_{domain}_{count}_sequences.fasta".format(
                    output=output, domain=domain_ID, count=len(arch_dict[domain_ID].keys())), 'a') as output_fasta:
                for key, values in arch_dict[domain_ID].items():
                    # print(values["sequence"])
                    output_fasta.write(">{key}\n{sequence}\n".format(key=key, sequence=values["sequence"]))


if __name__ == "__main__":
    list_of_domains, arch_dict = [], {}
    print("***** Input files parsing *****")
    text_file_parsing(args.domains, list_of_domains)
    domain_arch_parsing(args.arch, arch_dict)
    print("***** Domain sequences extraction *****")
    sequences_extraction(args.fasta, arch_dict)
    print("***** Output files writing *****")
    output_writing(args.output, arch_dict, list_of_domains)