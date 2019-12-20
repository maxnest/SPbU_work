try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_sign', type=argparse.FileType('r'), required=True)
args = parser.parse_args()


def nucl_parsing(nucl_list, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    nucl_list.extend([seq.id.split(" ")[0] for seq in contigs_nucl])


def amino_parsing(amino_list, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    amino_list.extend([seq.id.split(" ")[0].split(".p")[0] for seq in contigs_amino])


def significant(significant_list, table):
    for line in table:
        contig = line.strip().split(",")[0]
        significant_list.append(contig)


def amino_and_sign_writing(amino_list, significant_list, out, tag):
    all_ids = [contig for contig in significant_list]
    all_ids.extend(amino_list)
    print("In total, {count} sequences are assigned to '{tag}' category".format(count=len(set(all_ids)), tag=tag))
    with open("{out}.{tag}.txt".format(out=out, tag=tag), 'a') as output:
        for id in set(all_ids):
            output.write("{id}\n".format(id=id))


def without_protein_and_sign_writing(nucl_list, amino_list, significant_list, out, tag):
    all_ids = [contig for contig in significant_list]
    without_protein = []
    for contig in nucl_list:
        if contig not in amino_list:
            all_ids.append(contig)
            without_protein.append(contig)
    print("In total, {count} sequences are assigned to '{tag}' category, "
          "among them - {without} without predicted proteins".format(count=len(set(all_ids)), tag=tag,
                                                                     without=len(set(without_protein))))
    with open("{out}.{tag}.txt".format(out=out, tag=tag), 'a') as output:
        for id in set(all_ids):
            output.write("{id}\n".format(id=id))


if __name__ == "__main__":
    nucl_list, amino_list = [], []
    head_significant, tail_significant = [], []
    print("***** Input files parsing *****")
    nucl_parsing(nucl_list, args.nucl)
    amino_parsing(amino_list, args.amino)
    significant(head_significant, args.head_sign)
    significant(tail_significant, args.tail_sign)
    print("***** Output files writing *****")
    amino_and_sign_writing(amino_list, head_significant, "Pdum_ref_super", "head_significant_and_all_with_proteins")
    amino_and_sign_writing(amino_list, tail_significant, "Pdum_ref_super", "tail_significant_and_all_with_proteins")
    without_protein_and_sign_writing(nucl_list, amino_list, head_significant, "Pdum_ref_super",
                                     "head_significant_and_all_without_proteins")
    without_protein_and_sign_writing(nucl_list, amino_list, tail_significant, "Pdum_ref_super",
                                     "tail_significant_and_all_without_proteins")