try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--after_bam', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()


def write_fastq(out, id, seq, comment, qual, tag):
    if tag == "/1":
        output_name = "{out}.R1.fastq".format(out=out)
    else:
        output_name = "{out}.R2.fastq".format(out=out)

    with open(output_name, 'a') as output:
        output.write("{id}\n{seq}\n{comment}\n{qual}\n".format(id="{id}{tag}".format(id=id, tag=tag),
                                                               seq=seq, comment=comment, qual=qual))


for line in args.after_bam:
    R1_id = line.strip()
    R1_seq = next(args.after_bam).strip()
    R1_comment = next(args.after_bam).strip()
    R1_qual = next(args.after_bam).strip()
    write_fastq(args.out, R1_id, R1_seq, R1_comment, R1_qual, "/1")
    R2_id = next(args.after_bam).strip()
    R2_seq = next(args.after_bam).strip()
    R2_comment = next(args.after_bam).strip()
    R2_qual = next(args.after_bam).strip()
    write_fastq(args.out, R2_id, R2_seq, R2_comment, R2_qual, "/2")




