try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--sif', type=argparse.FileType('r'), required=True,
                    help="The simple interaction file (sif) with 3 columns: First_node Interaction Second_node")
parser.add_argument('--percent', type=str, required=True,
                    help="The percentage of top sequences to be included in the list of HUBs")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def graph_reconstruction(sif, graph_dict):
    header = sif.readline()
    for line in sif:
        description = line.strip().split("\t")
        first_node, second_node = description[0], description[2]

        if first_node not in graph_dict.keys():
            graph_dict[first_node] = {"in": [], "out": [], "description": []}

        if second_node not in graph_dict.keys():
            graph_dict[second_node] = {"in": [], "out": [], "description": []}

        if first_node not in graph_dict[second_node]["out"] and second_node not in graph_dict[first_node]["in"]:
            graph_dict[first_node]["out"].append(second_node)
            graph_dict[second_node]["in"].append(first_node)


def graph_analysis(graph_dict, percent):
    total_edges_counts, top_counts = [], []
    hubs, periphery, intermediate_nodes = 0, 0, 0
    for node, values in graph_dict.items():
        if len(set(values["in"] + values["out"])) not in total_edges_counts:
            total_edges_counts.append(len(set(values["in"] + values["out"])))

    top_indexes = sorted(range(len(total_edges_counts)), key=lambda i: total_edges_counts[i],
                         reverse=True)[:round((len(total_edges_counts)/100)*int(percent))]
    for index in top_indexes:
        top_counts.append(total_edges_counts[index])

    for node, values in graph_dict.items():
        if len(set(values["in"])) != 0 and len(set(values["out"])) == 0:
            values["description"].append("Periphery")
            periphery += 1
        elif len(set(values["in"] + values["out"])) in top_counts:
            values["description"].append("HUB")
            hubs += 1
        elif len(set(values["in"])) != 0 and len(set(values["out"])) != 0 and \
                len(set(values["in"] + values["out"])) not in top_counts:
            values["description"].append("Intermediate_node")
            intermediate_nodes += 1

    print("***** In graph {hubs} HUB nodes, {periphery} peripheral nodes, "
          "and {intermediate} intermediate nodes were founded *****".format(hubs=hubs, periphery=periphery,
                                                                            intermediate=intermediate_nodes))


def output_writing(output, graph_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Node_ID\tInput_edge_count\tOutput_edge_count\tTotal_edges_count\tDescription\n")
        for node, values in graph_dict.items():
            output_file.write("{id}\t{input}\t{output}\t{total}\t{description}\n".format(
                                id=node, input=len(set(values["in"])), output=len(set(values["out"])),
                                total=len(set(values["in"] + values["out"])),
                                description="|".join(values["description"])))


if __name__ == "__main__":
    graph_dict = {}
    graph_reconstruction(args.sif, graph_dict)
    graph_analysis(graph_dict, args.percent)
    output_writing(args.output, graph_dict)