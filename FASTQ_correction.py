try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fastq', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()

with open("{out}.fastq".format(out=args.out), 'a') as out_fastq:
    for line in args.fastq:
        name = line.strip()
        sequence = next(args.fastq).strip()
        next(args.fastq)
        quality = next(args.fastq).strip()
        if " /" in name:
            description = name.split(" /")
            new_name = "{id}/{num}".format(id=description[0][:-12], num=description[1])
            out_fastq.write("{name}\n{seq}\n+\n{score}\n".format(name=new_name, seq=sequence, score=quality))
        else:
            out_fastq.write("{name}\n{seq}\n+\n{score}\n".format(name=name, seq=sequence, score=quality))

