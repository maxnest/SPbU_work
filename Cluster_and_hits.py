try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--fd95_hits', type=argparse.FileType('r'), required=True)
parser.add_argument('--x4_hits', type=argparse.FileType('r'), required=True)
parser.add_argument('--x5_hits', type=argparse.FileType('r'), required=True)
parser.add_argument('--tag', type=str, required=True)
args = parser.parse_args()


def read_orthogroups(table, dict):
    head = table.readline()
    for line in table:
        description = line.strip().split("\t")
        group = description[0]
        if len(description) == 4:
            dict[group] = {"FD95": description[1].split(","), "X4": description[2].split(","),
                           "X5": description[3].split(",")}
        elif len(description) == 3:
            dict[group] = {"FD95": description[1].split(","), "X4": description[2].split(","),
                           "X5": []}
        else:
            dict[group] = {"FD95": description[1].split(","), "X4": [], "X5": []}

    for group, values in dict.items():
        for key, value in values.items():
            if len(value) == 1:
                for el in value:
                    if el == '':
                        value.remove(el)


def read_hits(hits, dict, line_tag):
    head = hits.readline()
    for line in hits:
        description = line.strip().split("\t")
        gene, hit = description[0], description[1]
        if line_tag == "X4":
            dict[gene] = hit
        else:
            dict["{line}|{gene}".format(line=line_tag, gene=gene)] = hit


def parsing_hits_dict(all_hits, prot_list, hits_dict):
    for el in prot_list:
        if el in hits_dict.keys():
            all_hits.append(hits_dict[el])


def cluster_analysis(ortho_dict, fd95_hits, x4_hits, x5_hits, cluster_with_hits):
        for cluster, values in ortho_dict.items():
            all_hits = []
            parsing_hits_dict(all_hits, values["FD95"], fd95_hits)
            parsing_hits_dict(all_hits, values["X4"], x4_hits)
            parsing_hits_dict(all_hits, values["X5"], x5_hits)
            #print(all_hits)
            if len(all_hits) == 0:
                cluster_with_hits[cluster] = "no_hit"
            else:
                cluster_with_hits[cluster] = ",".join(set(all_hits))
                #cluster_with_hits[cluster] = set(all_hits)


def summary(ortho_dict, cluster_with_hits, summary):
    for cluster, hit in cluster_with_hits.items():
        if hit not in summary.keys():
                summary[hit] = {"All_clusters": [cluster], "All_3_genome": [],
                                "Only_FD95_and_X5": [], "Only_FD95_and_X4": [], "Only_X5_and_X4": [],
                                "Only_FD95": [], "Only_X4": [], "Only_X5": []}
        elif hit in summary.keys():
                summary[hit]["All_clusters"].append(cluster)

    for hit, values in summary.items():
        for cluster in values["All_clusters"]:
            if len(ortho_dict[cluster]["FD95"]) != 0 and len(ortho_dict[cluster]["X4"]) != 0 and len(ortho_dict[cluster]["X5"]) != 0:
                values["All_3_genome"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) != 0 and len(ortho_dict[cluster]["X4"]) != 0 and len(ortho_dict[cluster]["X5"]) == 0:
                values["Only_FD95_and_X4"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) != 0 and len(ortho_dict[cluster]["X4"]) == 0 and len(ortho_dict[cluster]["X5"]) != 0:
                values["Only_FD95_and_X5"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) == 0 and len(ortho_dict[cluster]["X4"]) != 0 and len(ortho_dict[cluster]["X5"]) != 0:
                values["Only_X5_and_X4"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) != 0 and len(ortho_dict[cluster]["X4"]) == 0 and len(ortho_dict[cluster]["X5"]) == 0:
                values["Only_FD95"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) == 0 and len(ortho_dict[cluster]["X4"]) != 0 and len(ortho_dict[cluster]["X5"]) == 0:
                values["Only_X4"].append(cluster)
            elif len(ortho_dict[cluster]["FD95"]) == 0 and len(ortho_dict[cluster]["X4"]) == 0 and len(ortho_dict[cluster]["X5"]) != 0:
                values["Only_X5"].append(cluster)


if __name__ == "__main__":
    ortho_dict, fd95_hits, x4_hits, x5_hits, cluster_with_hits = {}, {}, {}, {}, {}
    summary_dict = {}
    read_orthogroups(args.tab, ortho_dict)
    read_hits(args.fd95_hits, fd95_hits, 'FD95')
    read_hits(args.x4_hits, x4_hits, 'X4')
    read_hits(args.x5_hits, x5_hits, 'X5')
    cluster_analysis(ortho_dict, fd95_hits, x4_hits, x5_hits, cluster_with_hits)
    summary(ortho_dict, cluster_with_hits, summary_dict)
    with open("{tag}_summary.tsv".format(tag=args.tag), 'a') as output:
        output.write("Hit(s)\t#_of_clusters\t#_with_all_3_genome\t#_with_only_FD95_and_X4\t#_with_only_FD95_and_X5"
                     "\t#_with_only_X5_and_X4\t#_with_only_FD95\t#_with_only_X4\t#_with_only_X5\n")
        for hit, values in summary_dict.items():
            output.write("{hit}\t{total}\t{all_3}\t{fd95_x4}\t{fd95_x5}\t{x5_x4}\t{fd95}\t{x4}\t{x5}\n".format(
                hit=hit, total=len(values["All_clusters"]), all_3=len(values["All_3_genome"]),
                fd95_x4=len(values["Only_FD95_and_X4"]), fd95_x5=len(values["Only_FD95_and_X5"]),
                x5_x4=len(values["Only_X5_and_X4"]), fd95=len(values["Only_FD95"]), x4=len(values["Only_X4"]),
                x5=len(values["Only_X5"]))
            )


    #with open("{tag}.tsv".format(tag=args.tag), 'a') as output:
    #    for key, value in cluster_with_hits.items():
    #        output.write("{cluster}\t{hits}\n".format(cluster=key, hits=",".join(value) if value != "no_hit" else "no_hit"))

