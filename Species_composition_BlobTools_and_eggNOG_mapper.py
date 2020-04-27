try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--good_fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta file with 'good' sequences (after TransRate v101)")
parser.add_argument('--blob_uniprot_tab', type=argparse.FileType('r'), required=True,
                    help="BlobTools table based on comparison with UniProt")
parser.add_argument('--blob_uniprot_hits', type=argparse.FileType('r'), required=True,
                    help="Table with unique taxon names and groups (For instance: Acidobacteria\tBacteria)")
#   parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True,
#                    help="Table with annotation results")
parser.add_argument('--wanted', type=argparse.FileType('r'), required=True,
                    help="Text file with IDs of wanted groups (one per line) (For instance: Metazoa)")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def fasta_parsing(good_fasta, good_fasta_dict):
    contigs = SeqIO.parse(good_fasta, 'fasta')
    for seq in contigs:
        good_fasta_dict[seq.id.split(" ")[0]] = {"Blob_UniProt": [], "eggNOG": []}


def blob_parsing(table, hits, tag, good_fasta_dict):
    hits_dict = {}
    for line in hits:
        description = line.strip().split("\t")
        hits_dict[description[0]] = description[1]

    for line in table:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            contig_ID, phylum = description[0], description[5]
            if phylum in hits_dict.keys():
                good_fasta_dict[contig_ID][tag].append(hits_dict[phylum])


#   def eggnog_parsing(table, good_fasta_dict):
#       for line in table:
#           if not line.startswith("#"):
#               description = line.strip().split("\t")
#               contig, annotation = description[0].split(".p")[0], description[1:]
#               if len(annotation) >= 16 and len(annotation[16]) != 0:
#                    good_fasta_dict[contig]["eggNOG"].append(annotation[16])


def wanted_contigs(good_fasta_dict, wanted, output):
    wanted_groups = []
    for line in wanted:
        wanted_groups.append(line.strip())

    for wanted_group in wanted_groups:
        with open("{output}.contigs_from_{wanted}.txt".format(output=output, wanted=wanted_group), 'a') as output_file:
            for contig, values in good_fasta_dict.items():
                if wanted_group == values["Blob_UniProt"][0]:
                    output_file.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    good_fasta_dict = {}
    fasta_parsing(args.good_fasta, good_fasta_dict)
    blob_parsing(args.blob_uniprot_tab, args.blob_uniprot_hits, "Blob_UniProt", good_fasta_dict)
#    eggnog_parsing(args.eggnog, good_fasta_dict)
    for contig, values in good_fasta_dict.items():

        if len(values["Blob_UniProt"]) == 0:
            values["Blob_UniProt"].append("-")

#       if len(values["eggNOG"]) == 0:
#            values["eggNOG"].append("-")
#       print("{contig}\t{values}".format(contig=contig, values=values))
    wanted_contigs(good_fasta_dict, args.wanted, args.output)
    