try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--kegg', type=str, required=True, help="The ID of the analyzed pathway in the KEGG.\n"
                                                            "For instance: ko04310 or ko04350")
parser.add_argument('--out', type=str, required=True, help="Prefix for output files")
args = parser.parse_args()


def table_parsing(tab, contig_dict):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        contig, pathway, sites_significant, sites_head_cluster, sites_tail_cluster, head_significant, head_cluster, \
        tail_significant, tail_cluster = description[0],  description[12], description[23], description[24], \
                                         description[25], description[26], description[27], description[28], \
                                         description[29]
        contig_dict[contig] = {"pathway": pathway, "sites_significant": sites_significant,
                               "sites_head_cluster": sites_head_cluster, "sites_tail_cluster": sites_tail_cluster,
                               "head_significant": head_significant, "head_cluster": head_cluster,
                               "tail_significant": tail_significant, "tail_cluster": tail_cluster}


def append_contig_to_significant(kegg_dict, contig, values, values_tag, kegg_tag):
    if values[values_tag] != '-':
        kegg_dict[kegg_tag].append(contig)


def append_contig_to_cluster(kegg_dict, contig, values, values_tag):
    if values[values_tag] != '-':
        if values[values_tag] not in kegg_dict.keys():
            kegg_dict[values[values_tag]] = []
        kegg_dict[values[values_tag]].append(contig)


def kegg_summary(contig_dict, kegg, kegg_dict):
    for contig, values in contig_dict.items():
        if kegg in values["pathway"].split(","):
            kegg_dict["total"].append(contig)
            append_contig_to_significant(kegg_dict, contig, values, "sites_significant", "sites")
            append_contig_to_significant(kegg_dict, contig, values, "head_significant", "head")
            append_contig_to_significant(kegg_dict, contig, values, "tail_significant", "tail")
            append_contig_to_cluster(kegg_dict, contig, values, "sites_head_cluster")
            append_contig_to_cluster(kegg_dict, contig, values, "sites_tail_cluster")
            append_contig_to_cluster(kegg_dict, contig, values, "head_cluster")
            append_contig_to_cluster(kegg_dict, contig, values, "tail_cluster")


def output_writing(out, kegg, kegg_dict):
    cluster_keys = [cluster for cluster in kegg_dict.keys() if cluster not in ["total", "sites", "head", "tail"]]

    with open("{out}.{kegg}_summary.tsv".format(out=out, kegg=kegg), 'a') as output:
        output.write("### In total, {count} sequences are assigned to {kegg} pathway\n"
                     "### Among them:\n"
                     "### {sites} were previously classified as sites-significant:\t{sites_contigs}\n"
                     "### {head} were previously classified as head-significant:\t{head_contigs}\n"
                     "### {tail} were previously classified as tail-significant:\t{tail_contigs}\n"
                     "Cluster\tNumber of 'pathway'-contigs in cluster\tContigs (comma separated)\n".format(
                      count=len(set(kegg_dict["total"])), kegg=kegg,
                      sites=len(set(kegg_dict["sites"])), sites_contigs=",".join(set(kegg_dict["sites"])),
                      head=len(set(kegg_dict["head"])), head_contigs=",".join(set(kegg_dict["head"])),
                      tail=len(set(kegg_dict["tail"])), tail_contigs=",".join(set(kegg_dict["tail"]))
        ))
        for cluster in cluster_keys:
            output.write("{cluster}\t{length}\t{contigs}\n".format(cluster=cluster, length=len(set(kegg_dict[cluster])),
                                                                   contigs=",".join(set(kegg_dict[cluster]))))


if __name__ == "__main__":
    contig_dict, kegg_dict = {}, {"total": [], "sites": [], "head": [], "tail": []}
    print("***** Input table parsing *****")
    table_parsing(args.tab, contig_dict)
    print("***** Search for sequences related to {kegg} *****".format(kegg=args.kegg))
    kegg_summary(contig_dict, args.kegg, kegg_dict)
    print("***** Output file writing *****")
    output_writing(args.out, args.kegg, kegg_dict)