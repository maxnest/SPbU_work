try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True, help="Pdum_summary.py output file")
parser.add_argument('--output', type=str, required=True, help="Prefix for output files")
args = parser.parse_args()


def read_tab(tab, contig_dict):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        contig, nt, swiss, go_bio, go_mol, head_cluster, tail_cluster = description[0], description[1], \
                      description[3], description[5], description[6], description[11], description[13]
        contig_dict[contig] = {"nt": nt, "swiss": swiss, "go_bio": go_bio.split("|")[0], "go_mol": go_mol.split("|")[0],
                               "head_cluster": head_cluster, "tail_cluster": tail_cluster}


def clusters_grouping(contig_dict, cluster_dict, tag):
    for contig, values in contig_dict.items():
        if values["{tag}_cluster".format(tag=tag)] != "-":
            cluster_ID = values["{tag}_cluster".format(tag=tag)]
            if values["{tag}_cluster".format(tag=tag)] not in cluster_dict.keys():
                cluster_dict[cluster_ID] = {"contigs": [], "nt_hits": [], "swiss_hits": [], "go_bio": [], "go_mol": []}
                cluster_dict[cluster_ID]["contigs"].append(contig)
                cluster_dict[cluster_ID]["nt_hits"].append(values["nt"])
                cluster_dict[cluster_ID]["swiss_hits"].append(values["swiss"])
                cluster_dict[cluster_ID]["go_bio"].append(values["go_bio"])
                cluster_dict[cluster_ID]["go_mol"].append(values["go_mol"])
            else:
                cluster_dict[cluster_ID]["contigs"].append(contig)
                cluster_dict[cluster_ID]["nt_hits"].append(values["nt"])
                cluster_dict[cluster_ID]["swiss_hits"].append(values["swiss"])
                cluster_dict[cluster_ID]["go_bio"].append(values["go_bio"])
                cluster_dict[cluster_ID]["go_mol"].append(values["go_mol"])


def summary(cluster_dict, tag, output):
    with open("{output}.{tag}_clusters.nohits.tsv".format(output=output, tag=tag), 'a') as output_file:
        output_file.write("Cluster_ID\tCluster_size\tNCBInt (no-hits %)\tSwiss (without_ORFS %|no-hits %)\n")
        for cluster, values in cluster_dict.items():
            output_file.write("{id}\t{size}\t{ncbi_nohit}\t{swiss_without_orf}|{swiss_nohit}\n".format(
                id=cluster, size=len(values["contigs"]),
                ncbi_nohit=round((values["nt_hits"].count("No hit")/len(values["nt_hits"]))*100),
                swiss_without_orf=round((values["swiss_hits"].count("-")/len(values["swiss_hits"]))*100),
                swiss_nohit=round((values["swiss_hits"].count("No hit")/len(values["swiss_hits"]))*100)
            ))

    for cluster, values in cluster_dict.items():
        with open("{output}.{cluster}.GO_bio_counts.tsv".format(output=output, cluster=cluster), 'a') as go_bio_counts:
            for go_bio in set(values["go_bio"]):
                if go_bio != "-":
                    go_bio_counts.write("{ID}\t{counts}\n".format(ID=go_bio, counts=values["go_bio"].count(go_bio)))

    for cluster, values in cluster_dict.items():
        with open("{output}.{cluster}.GO_mol_counts.tsv".format(output=output, cluster=cluster, tag=tag), 'a') \
                as go_mol_counts:
            for go_mol in set(values["go_mol"]):
                if go_mol != "-":
                    go_mol_counts.write("{ID}\t{counts}\n".format(ID=go_mol, counts=values["go_mol"].count(go_mol)))


if __name__ == "__main__":
    contig_dict, head_clusters, tail_clusters = {}, {}, {}
    print("***** Input file parsing *****")
    read_tab(args.tab, contig_dict)
    print("***** Clusters description parsing ******")
    clusters_grouping(contig_dict, head_clusters, "head")
    clusters_grouping(contig_dict, tail_clusters, "tail")
    print("***** Output files creating *****")
    summary(head_clusters, "head", args.output)
    summary(tail_clusters, "tail", args.output)