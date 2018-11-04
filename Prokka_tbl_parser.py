try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--prokka_tbl', type=argparse.FileType('r'), required=True)
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str)
args = parser.parse_args()
fasta_lib = SeqIO.parse(args.fasta, "fasta")


prokka_dict, fasta_dict = {}, {}
contigs_keys = []

contigs_with_hits = []
contigs_with_predicted_only = []
contigs_without_hits = []


def write_fasta(fasta_dict, list, output, pattern):
    with open("{output}_{pattern}.fasta".format(output=output, pattern=pattern), 'a') as output_file:
        for el in list:
            output_file.write("{name}\n{seq}\n".format(name=">{el}".format(el=el), seq=fasta_dict[el]))


for line in args.prokka_tbl:
    description = line.strip().split("\t")
    if description[0].startswith(">"):
        contig_id = line.strip().split(" ")[1]
        prokka_dict[contig_id] = []
        contigs_keys.append(contig_id)
    else:
        if "locus_tag" in description:
            product = next(args.prokka_tbl).strip().split("\t")[-1]
            prokka_dict[contigs_keys[-1]].append({description[-1]: product})


for contig_id in prokka_dict.keys():
    if len(prokka_dict[contig_id]) == 0:
        contigs_without_hits.append(contig_id)
    else:
        products = []
        for locus in prokka_dict[contig_id]:
            for product in locus.values():
                products.append(product)
        if len(set(products)) == 1 and "hypothetical protein" in set(products):
            contigs_with_predicted_only.append(contig_id)
        else:
            contigs_with_hits.append(contig_id)

print("# of contigs with hits: {len_hits}\n"
      "# of contigs with predicted only: {len_predicted}\n"
      "# of contigs without hits: {len_without}\n".format(len_hits=len(contigs_with_hits),
                                                          len_predicted=len(contigs_with_predicted_only),
                                                          len_without=len(contigs_without_hits)))

for fasta in fasta_lib:
    name, sequence = fasta.id, fasta.seq.tostring()
    fasta_dict[name] = sequence

write_fasta(fasta_dict, contigs_with_hits, args.output, "with_Prokka_hits")
write_fasta(fasta_dict, contigs_without_hits, args.output, "without_Prokka_hits")
write_fasta(fasta_dict, contigs_with_predicted_only, args.output, "with_Prokka_predicted_only")



