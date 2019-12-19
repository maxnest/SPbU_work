try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True, help="Prefix for output files")
args = parser.parse_args()


def table_parsing(tab, contig_dict):
    header = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        contig, nt, nr, swiss, domains, go, ec, ko, pathway, module, reaction, rclass, brite, tc, cazy, bigg, og, cog, \
        sites_significant, sites_head_cluster, sites_tail_cluster, head_significant, head_cluster, tail_significant, \
        tail_cluster = description[0], description[2], description[4], description[6], description[8], description[9], \
                       description[10], description[11], description[12], description[13], description[14], \
                       description[15], description[16], description[17], description[18], description[19], \
                       description[20], description[21], description[23], description[24], description[25], \
                       description[26], description[27], description[28], description[29]

        contig_dict[contig] = {"nt": nt, "nr": nr, "swiss": swiss, "domains": domains, "go": go, "ec": ec, "ko": ko,
                               "pathway": pathway, "module": module, "reaction": reaction, "rclass": rclass,
                               "brite": brite, "tc": tc, "cazy": cazy, "bigg": bigg, "og": og, "cog": cog,
                               "sites_significant": sites_significant, "sites_head_cluster": sites_head_cluster,
                               "sites_tail_cluster": sites_tail_cluster, "head_significant": head_significant,
                               "head_cluster": head_cluster, "tail_significant": tail_significant,
                               "tail_cluster": tail_cluster}


def info_appending_for_significant(contig, values, tag, dict, keys):
    if values[tag] != '-':
        dict["all"].append(contig)
        for key in keys:
            if values[key] != "-" and values[key] != "No hit":
                dict["annotated"].append(contig)


def info_appending_for_clusters(contig, values, tag, dict, keys):
    if values[tag] != '-':
        if values[tag] not in dict.keys():
            dict[values[tag]] = {"all": [], "annotated": [], "cogs": []}
        dict[values[tag]]["all"].append(contig)
        for key in keys:
            if values[key] != '-' and values[key] != "No hit":
                dict[values[tag]]["annotated"].append(contig)

            if key == "cog" and values[key] != "-":
                dict[values[tag]]["cogs"].extend(list(values[key]))


def annotation_summary(contig_dict, sites_significant, head_significant, tail_significant,
                       sites_head_clusters, sites_tail_clusters, head_clusters, tail_clusters):
    keys = ["nt", "nr", "swiss", "domains", "go", "ec", "ko", "pathway", "module", "reaction", "rclass",
            "brite", "tc", "cazy", "bigg", "og", "cog"]
    for contig, values in contig_dict.items():
        info_appending_for_significant(contig, values, "sites_significant", sites_significant, keys)
        info_appending_for_significant(contig, values, "head_significant", head_significant, keys)
        info_appending_for_significant(contig, values, "tail_significant", tail_significant, keys)
        info_appending_for_clusters(contig, values, "sites_head_cluster", sites_head_clusters, keys)
        info_appending_for_clusters(contig, values, "sites_tail_cluster", sites_tail_clusters, keys)
        info_appending_for_clusters(contig, values, "head_cluster", head_clusters, keys)
        info_appending_for_clusters(contig, values, "tail_cluster", tail_clusters, keys)
        for key in keys:
            if values[key] != '-' and values[key] != "No hit":
                annotated_contigs.append(contig)


def annotation_summary_output_writing(contig_dict, annotated_contigs, sites_significant, head_significant,
                                      tail_significant, sites_head_clusters, sites_tail_clusters, head_clusters,
                                      tail_clusters, out):
    with open("{out}.annotation_results_summary.txt".format(out=out), 'a') as output:
        output.write("Number of annotated contigs: {contigs} ({contigs_percent}%)\n"
                     "Number of annotated sites_significant: {sites} ({sites_percent}%)\n"
                     "Number of annotated head_significant: {head} ({head_percent}%)\n"
                     "Number of annotated tail_significant: {tail} ({tail_percent}%)\n"
                     "Sites_head_clusters:\n{sites_head_clusters}\n"
                     "Sites_tail_clusters:\n{sites_tail_clusters}\n"
                     "Head_clusters:\n{head_clusters}\n"
                     "Tail_clusters:\n{tail_clusters}\n".format(
                      contigs=len(set(annotated_contigs)),
                      contigs_percent=(len(set(annotated_contigs))/len(contig_dict.keys()))*100,
                      sites=len(set(sites_significant["annotated"])),
                      sites_percent=(len(set(sites_significant["annotated"]))/len(set(sites_significant["all"])))*100,
                      head=len(set(head_significant["annotated"])),
                      head_percent=(len(set(head_significant["annotated"]))/len(set(head_significant["all"])))*100,
                      tail=len(set(tail_significant["annotated"])),
                      tail_percent=(len(set(tail_significant["annotated"]))/len(set(tail_significant["all"])))*100,
                      sites_head_clusters="\n".join(["{cluster}:{annotated}({percent}%)".format(
                          cluster=cluster, annotated=len(set(sites_head_clusters[cluster]["annotated"])),
                          percent=(len(set(sites_head_clusters[cluster]["annotated"]))/
                                   len(set(sites_head_clusters[cluster]["all"])))*100)
                          for cluster in sites_head_clusters.keys()]),
                      sites_tail_clusters="\n".join(["{cluster}:{annotated}({percent}%)".format(
                          cluster=cluster, annotated=len(set(sites_tail_clusters[cluster]["annotated"])),
                          percent=(len(set(sites_tail_clusters[cluster]["annotated"]))/
                                   len(set(sites_tail_clusters[cluster]["all"])))*100)
                          for cluster in sites_tail_clusters.keys()]),
                      head_clusters="\n".join(["{cluster}:{annotated}({percent}%)".format(
                          cluster=cluster, annotated=len(set(head_clusters[cluster]["annotated"])),
                          percent=(len(set(head_clusters[cluster]["annotated"]))/
                                   len(set(head_clusters[cluster]["all"])))*100)
                          for cluster in head_clusters.keys()]),
                      tail_clusters="\n".join(["{cluster}:{annotated}({percent}%)".format(
                          cluster=cluster, annotated=len(set(tail_clusters[cluster]["annotated"])),
                          percent=(len(set(tail_clusters[cluster]["annotated"]))/
                                   len(set(tail_clusters[cluster]["all"])))*100)
                          for cluster in tail_clusters.keys()])))


def cluster_annotation_summary(sites_head_clusters, sites_tail_clusters, head_clusters, tail_clusters, out):
    COG_keys = ["D:Cell cycle control, cell division, chromosome partitioning",
                "M:Cell wall/membrane/envelope biogenesis",
                "N:Cell motility",
                "O:Post-translational modification, protein turnover, and chaperones",
                "T:Signal transduction mechanisms",
                "U:Intracellular trafficking, secretion, and vesicular transport",
                "V:Defense mechanisms",
                "W:Extracellular structures",
                "Y:Nuclear structure",
                "Z:Cytoskeleton",
                "A:RNA processing and modification",
                "B:Chromatin structure and dynamics",
                "J:Translation, ribosomal structure and biogenesis",
                "K:Transcription",
                "L:Replication, recombination and repair",
                "C:Energy production and conversion",
                "E:Amino acid transport and metabolism",
                "F:Nucleotide transport and metabolism",
                "G:Carbohydrate transport and metabolism",
                "H:Coenzyme transport and metabolism",
                "I:Lipid transport and metabolism",
                "P:Inorganic ion transport and metabolism",
                "Q:Secondary metabolites biosynthesis, transport, and catabolism",
                "R:General function prediction only",
                "S:Function unknown"]
    sites_head_keys, sites_tail_keys, head_keys, tail_keys = [key for key in sites_head_clusters.keys()], \
                                                             [key for key in sites_tail_clusters.keys()], \
                                                             [key for key in head_clusters.keys()], \
                                                             [key for key in tail_clusters.keys()]
    with open("{out}.cluster_COG_annotation_summary.tsv".format(out=out), 'a') as output:
        output.write("COG\tDescription\t{head}\t{tail}\t{sites_head}\t{sites_tail}\n".format(head="\t".join(head_keys),
               tail="\t".join(tail_keys), sites_head="\t".join(sites_head_keys), sites_tail="\t".join(sites_tail_keys)))
        for COG in COG_keys:
            output.write("{COG}\t{Description}\t{head}\t{tail}\t{sites_head}\t{sites_tail}\n".format(
                COG=COG.split(":")[0], Description=COG.split(":")[1],
                head="\t".join([str(head_clusters[key]["cogs"].count(COG.split(":")[0])) for key in head_keys]),
                tail="\t".join([str(tail_clusters[key]["cogs"].count(COG.split(":")[0])) for key in tail_keys]),
                sites_head="\t".join([str(sites_head_clusters[key]["cogs"].count(COG.split(":")[0]))
                                      for key in sites_head_keys]),
                sites_tail="\t".join([str(sites_tail_clusters[key]["cogs"].count(COG.split(":")[0]))
                                      for key in sites_tail_keys])
            ))


if __name__ == "__main__":
    contig_dict = {}
    annotated_contigs, sites_significant, head_significant, tail_significant, = [], \
    {"all": [], "annotated": []}, {"all": [], "annotated": []}, {"all": [], "annotated": []}
    sites_head_clusters, sites_tail_clusters, head_clusters, tail_clusters = {}, {}, {}, {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, contig_dict)
    print("***** Annotation summary *****")
    annotation_summary(contig_dict, sites_significant, head_significant, tail_significant, sites_head_clusters,
                       sites_tail_clusters, head_clusters, tail_clusters)
    annotation_summary_output_writing(contig_dict, annotated_contigs, sites_significant, head_significant,
                                      tail_significant, sites_head_clusters, sites_tail_clusters, head_clusters,
                                      tail_clusters, args.out)
    print("***** Summary of COG-annotation of clusters *****")
    cluster_annotation_summary(sites_head_clusters, sites_tail_clusters, head_clusters, tail_clusters, args.out)