try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--blastx_ref', type=argparse.FileType('r'), required=True,
                    help="BLASTx table in 6 outfmt with results of comparison contigs with reference proteome")
parser.add_argument('--blastx_uniprot', type=argparse.FileType('r'), required=True,
                    help="BLASTx table in 6 outfmt with results of comparison contigs with UniProt database")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def blast_tab_parsing(blast_tab, blast_dict):
    for line in blast_tab:
        description = line.strip().split("\t")
        qseqid, sseqid, pident, length = description[0], description[1], float(description[2]), float(description[3])
        if qseqid not in blast_dict:
            blast_dict[qseqid] = []

        blast_dict[qseqid].append({"sseqid": sseqid, "pident": pident, "length": length})


def contigs_selection(contigs_dict, ref_dict, uniprot_dict, reliable):
    for contig in contigs_dict.keys():
        if contig in ref_dict.keys() and contig in uniprot_dict.keys():
            ref_pident = [el["pident"] for el in ref_dict[contig]]
            uniprot_pident = [el["pident"] for el in uniprot_dict[contig]]
            ref_length = [el["length"] for el in ref_dict[contig]]
            uniprot_length = [el["length"] for el in uniprot_dict[contig]]
            if max(ref_pident) > max(uniprot_pident) or max(ref_length) > max(uniprot_length):
                reliable.append(contig)
    print("### Number of reliable contigs: {len} ###".format(len=len(set(reliable))))


def output_writing(output, reliable, ref_dict, uniprot_dict):
    with open("{output}.IDs_of_all_reliable_contigs.txt".format(output=output), 'a') as reliable_txt:
        for contig in set(reliable):
            reliable_txt.write("{contig}\n".format(contig=contig))

    with open("{output}.all_reliable_contigs_alignment_results.tsv".format(output=output), 'a') as reliable_table:
        reliable_table.write("Contig_ID\tReference:max_pident\tReference:max_align_length\t"
                             "UniProt:max_pident\tUniProt:max_align_length\n")
        for contig in set(reliable):
            ref_pident = [el["pident"] for el in ref_dict[contig]]
            uniprot_pident = [el["pident"] for el in uniprot_dict[contig]]
            ref_length = [el["length"] for el in ref_dict[contig]]
            uniprot_length = [el["length"] for el in uniprot_dict[contig]]

            reliable_table.write("{id}\t{ref_pident}\t{ref_len}\t{uniprot_pident}\t{uniprot_len}\n".format(
                id=contig, ref_pident=max(ref_pident), ref_len=max(ref_length),
                uniprot_pident=max(uniprot_pident), uniprot_len=max(uniprot_length)
            ))


if __name__ == "__main__":
    contigs_dict, reliable, ref_dict, uniprot_dict = {}, [], {}, {}

    contigs_fasta = SeqIO.parse(args.fasta, 'fasta')
    for seq in contigs_fasta:
        contigs_dict[seq.id] = seq.seq

    blast_tab_parsing(args.blastx_ref, ref_dict)
    blast_tab_parsing(args.blastx_uniprot, uniprot_dict)

    contigs_selection(contigs_dict, ref_dict, uniprot_dict, reliable)
    output_writing(args.output, reliable, ref_dict, uniprot_dict)