try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--r1', type=argparse.FileType('r'), required=True)
parser.add_argument('--r2', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def jaccard_similarity(one_list, other_list):
    """
    FROM: http://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
    The function measure the similarity between two sets of genes (Jaccard similarity index)
    @param one_list: list of genes
    @param other_list: list of genes
    """
    intersection_cardinality = len(set.intersection(*[set(one_list), set(other_list)]))
    union_cardinality = len(set.union(*[set(one_list), set(other_list)]))
    if union_cardinality != 0:
        return round(intersection_cardinality / float(union_cardinality), 2)
    else:
        return 0


def Jaccard_for_lib_pair(r1_fastq, r2_fastq, out):
    r1_ids, r2_ids = [read.id[:-1] for read in r1_fastq], [read.id[:-1] for read in r2_fastq]
    print("Lengths: IDs_of_R1 - {len1}; IDs_of_R2 - {len2}".format(len1=len(r1_ids), len2=len(r2_ids)))
    print("First lines: R1 - {fl_1}; R2 - {fl_2}".format(fl_1=r1_ids[0], fl_2=r2_ids[0]))
    print("Last lines: R1 - {ll_1}; R2 - {ll_2}".format(ll_1=r1_ids[-1], ll_2=r2_ids[-1]))
    print("Jaccard:{out} - {value}".format(out=out, value=jaccard_similarity(r1_ids, r2_ids)))


if __name__ == "__main__":
    r1_fastq, r2_fastq = SeqIO.parse(args.r1, "fastq"), SeqIO.parse(args.r2, "fastq")
    Jaccard_for_lib_pair(r1_fastq, r2_fastq, args.out)
