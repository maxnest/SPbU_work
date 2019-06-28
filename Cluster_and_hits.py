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
                cluster_with_hits[cluster] = set(all_hits)


if __name__ == "__main__":
    ortho_dict, fd95_hits, x4_hits, x5_hits, cluster_with_hits = {}, {}, {}, {}, {}
    read_orthogroups(args.tab, ortho_dict)
    read_hits(args.fd95_hits, fd95_hits, 'FD95')
    read_hits(args.x4_hits, x4_hits, 'X4')
    read_hits(args.x5_hits, x5_hits, 'X5')
    cluster_analysis(ortho_dict, fd95_hits, x4_hits, x5_hits, cluster_with_hits)

    with open("{tag}.tsv".format(tag=args.tag), 'a') as output:
        for key, value in cluster_with_hits.items():
            output.write("{cluster}\t{hits}\n".format(cluster=key, hits=",".join(value) if value != "no_hit" else "no_hit"))

