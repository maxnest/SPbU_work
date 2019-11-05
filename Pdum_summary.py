try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nt', type=argparse.FileType('r'), required=True)
parser.add_argument('--GO', type=argparse.FileType('r'), required=True)
parser.add_argument('--sites_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ID": [seq.id.split(" ")[0]], "nucl": seq.seq, "amino": "",
                               "domains": [], "GO_bio": [], "GO_mol": [], "GO_cell": [],
                               "Anno_swiss": {"hit": [], "identity": []}, "Anno_nt": {"hit": [], "identity": []},
                               "head_significant": [], "tail_significant": [], "sites_significant": [],
                               "head_cluster": [], "tail_cluster": []}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".p")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"


def GO_parsing(contig_dict, GO):
    head = GO.readline()
    for line in GO:
        description = line.strip().split("\t")
        contig_ID, category, ID, rank, anno = description[0].strip().split(".p")[0], description[1], \
                                              "GO:{ID}".format(ID=description[2]), description[-1], description[3]
        if contig_ID in contig_dict.keys():
            if rank == "1":
                if category == "MF":
                    contig_dict[contig_ID]["GO_mol"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
                elif category == "CC":
                    contig_dict[contig_ID]["GO_cell"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
                elif category == "BP":
                    contig_dict[contig_ID]["GO_bio"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))

    for contig, values in contig_dict.items():
        if len(values["GO_mol"]) == 0:
            contig_dict[contig]["GO_mol"].append("-")

        if len(values["GO_cell"]) == 0:
            contig_dict[contig]["GO_cell"].append("-")

        if len(values["GO_bio"]) == 0:
            contig_dict[contig]["GO_bio"].append("-")


def BLAST_annotation(contig_dict, BLAST, key_tag):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        if key_tag == "swiss":
            ID, hit_name = description[0].strip().split(".p")[0], description[3]
        else:
            ID, hit_name = description[0].strip().split(" ")[0], description[3]

        if ID in contig_dict.keys():
            if hit_name == "no hits":
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
            else:
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append(description[8])

    if key_tag == "swiss":
        for contig in contig_dict.keys():
            if contig_dict[contig]["amino"] == "without ORF":
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["hit"].append("-")
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["identity"].append("-")


def domains_parsing(contig_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        ID, domain_name = description[0].strip().split(".p")[0], description[1]
        contig_dict[ID]["domains"].append(domain_name)

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["domains"]) == 0:
            contig_dict[contig]["domains"].append("-")


def clusters(contig_dict, table, dict, tag):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)

    for contig in contig_dict.keys():
        for cluster, genes in dict.items():
            if contig in genes:
                contig_dict[contig]["{tag}_cluster".format(tag=tag)].append(cluster)

    for contig, values in contig_dict.items():
        if len(values["{tag}_cluster".format(tag=tag)]) == 0:
            values["{tag}_cluster".format(tag=tag)].append("-")


def significant(contig_dict, table, tag):
    for line in table:
        contig = line.strip().split(",")[0]
        contig_dict[contig]["{tag}_significant".format(tag=tag)].append("*")

    for contig, values in contig_dict.items():
        if len(values["{tag}_significant".format(tag=tag)]) == 0:
            values["{tag}_significant".format(tag=tag)].append("-")


def write_output_files(contig_dict, output):
    # NB! Without AA sequences now
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Contig_ID\tNCBInt\tNCBInt_identity\tSwiss\tSwiss_identity\tGO_bio\tGO_mol\tGO_cell\t"
                          "Domains_arch\tSites_significant\tHead_significant\tHead_cluster\tTail_significant\t"
                          "Tail_cluster\n")

        for contig, values in contig_dict.items():
            output_file.write("{id}\t{nt}\t{nt_identity}\t{swiss}\t{swiss_identity}\t{bio}\t{mol}\t{cell}\t"
                              "{domains}\t{sites_sign}\t{head_sign}\t{head_cluster}\t{tail_sign}\t{tail_cluster}\n".format(
                id=contig, nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                swiss=values["Anno_swiss"]["hit"][0], swiss_identity=values["Anno_swiss"]["identity"][0],
                bio=values["GO_bio"][0], mol=values["GO_mol"][0], cell=values["GO_cell"][0],
                domains="|".join(values["domains"]), sites_sign=values["sites_significant"][0],
                head_sign=values["head_significant"][0], head_cluster=values["head_cluster"][0],
                tail_sign=values["tail_significant"][0], tail_cluster=values["tail_cluster"][0]
            ))


if __name__ == "__main__":
    contig_dict = {}
    head_clusters, tail_clusters = {}, {}
    print("***** Fasta files parsing *****")
    nucl_parsing(contig_dict, args.nucl)
    amino_parsing(contig_dict, args.amino)
    print("***** GO parsing *****")
    GO_parsing(contig_dict, args.GO)
    print("***** BLAST results parsing *****")
    BLAST_annotation(contig_dict, args.blast_swiss, "swiss")
    BLAST_annotation(contig_dict, args.blast_nt, "nt")
    print("***** Domain architecture parsing *****")
    domains_parsing(contig_dict, args.domains)
    print("***** Cluster parsing *****")
    clusters(contig_dict, args.head_clusters, head_clusters, "head")
    clusters(contig_dict, args.tail_clusters, tail_clusters, "tail")
    print("***** Tables with IDs of significant sequences parsing *****")
    significant(contig_dict, args.sites_sign, "sites")
    significant(contig_dict, args.head_sign, "head")
    significant(contig_dict, args.tail_sign, "tail")
    print("***** Output file creating *****")
    write_output_files(contig_dict, args.output)

