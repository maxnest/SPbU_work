try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--r1', type=argparse.FileType('r'), required=True)
parser.add_argument('--r2', type=argparse.FileType('r'), required=True)
parser.add_argument('--bamfilter', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()


def fastq_parsing(fastq, fastq_dict):
    for line in fastq:
        ID, seq, comment, qual = line.strip().split("/")[0], next(fastq).strip(), next(fastq).strip(), next(fastq).strip()
        fastq_dict[ID] = {"seq": seq, "comment": comment, "qual": qual}


def bamfilter_parsing(bamfilter, IDs):
    for line in bamfilter:
        ID, seq, comment, qual = line.strip(), next(bamfilter).strip(), next(bamfilter).strip(), next(bamfilter).strip()
        IDs.append(ID)


def output_writing(out, IDs, r1_dict, r2_dict):
    Uniq_IDs = set(IDs)
    for ID in Uniq_IDs:
        with open("{out}.R1.fastq".format(out=out), 'a') as R1:
            R1.write("{id}\n{seq}\n{comment}\n{qual}\n".format(
                id="{id}/1".format(id=ID), seq=r1_dict[ID]["seq"], comment=r1_dict[ID]["comment"],
                qual=r1_dict[ID]["qual"]
            ))
        with open("{out}.R2.fastq".format(out=out), 'a') as R2:
            R2.write("{id}\n{seq}\n{comment}\n{qual}\n".format(
                id="{id}/2".format(id=ID), seq=r2_dict[ID]["seq"], comment=r2_dict[ID]["comment"],
                qual=r2_dict[ID]["qual"]
            ))


if __name__ == "__main__":
    r1_dict, r2_dict, IDs_list = {}, {}, []
    print("***** Input files parsing *****")
    fastq_parsing(args.r1, r1_dict)
    fastq_parsing(args.r2, r2_dict)
    bamfilter_parsing(args.bamfilter, IDs_list)
    print("***** Output files writing *****")
    output_writing(args.out, IDs_list, r1_dict, r2_dict)