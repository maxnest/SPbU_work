try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="SecretomeP output table with results and without any process description")
parser.add_argument('--out', type=str)
args = parser.parse_args()


def table_parsing(table, signalp, nonclass):
    for line in table:
        description = line.strip().split(" ")
        ID, nn, warning = description[0].split(".p")[0], 0, description[-1]

        for el in description[2:]:
            if len(el) != 0:
                nn += float(el.strip().split("\t")[0])
                break

        if warning == "SignalP":
            if ID not in signalp:
                if ID.endswith("."):
                    signalp.append(ID[:-1])
                else:
                    signalp.append(ID)
        else:
            if nn >= 0.9 and warning == "-":
                if ID not in nonclass:
                    if ID.endswith("."):
                        nonclass.append(ID[:-1])
                    else:
                        nonclass.append(ID)


def sorting(fasta_seqs, signalp, nonclass, signalp_dict, nonclass_dict):
    signalp_used, nonclass_used = [], []
    for fasta in fasta_seqs:
        name, sequence = fasta.id, fasta.seq

        if name.split("Cluster-")[1].split(".p")[0] in signalp:
            signalp_dict[name] = sequence
            signalp_used.append(name.split("Cluster-")[1].split(".p")[0])
        elif name.split("Cluster-")[1].split(".p")[0] in nonclass:
            nonclass_dict[name] = sequence
            nonclass_used.append(name.split("Cluster-")[1].split(".p")[0])


def output_writing(signalp_dict, nonclass_dict, out):
    with open("{out}.SignalP_from_SecretomeP.fasta".format(out=out), 'a') as signalp_fasta:
        for name, sequence in signalp_dict.items():
            signalp_fasta.write(">{name}\n{seq}\n".format(name=name, seq=sequence))

    with open("{out}.non-classical.fasta".format(out=out), 'a') as nc_fasta:
        for name, sequence in nonclass_dict.items():
            nc_fasta.write(">{name}\n{seq}\n".format(name=name, seq=sequence))


if __name__ == "__main__":
    fasta_seqs = SeqIO.parse(args.fasta, "fasta")
    signalp, nonclass = [], []
    signalp_dict, nonclass_dict = {}, {}
    print("***** SecretomeP output parsing *****")
    table_parsing(args.tab, signalp, nonclass)
    print("Length SignalP: {sp}; Non-classical: {nc}".format(sp=len(signalp), nc=len(nonclass) ))
    print("***** Sorting *****")
    sorting(fasta_seqs, signalp, nonclass, signalp_dict, nonclass_dict)
    print("***** Output files creating *****")
    output_writing(signalp_dict, nonclass_dict, args.out)
