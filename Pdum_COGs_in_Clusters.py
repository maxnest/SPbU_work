try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def table_parsing(tab, cluster_dict):
    for line in tab:
        if not line.startswith("#") and not line.startswith("Contig_ID"):
            description = line.strip().split("\t")
            contig_ID, COG, sites_head_cluster, sites_tail_cluster, head_cluster, tail_cluster = description[0], \
                            list(description[21]), description[24], description[25], description[27], description[29]
            clusters = [sites_head_cluster, sites_tail_cluster, head_cluster, tail_cluster]
            for cluster in clusters:
                if cluster != "-":
                    if cluster not in cluster_dict.keys():
                        cluster_dict[cluster] = {"contigs": [], "COGs": []}
                    cluster_dict[cluster]["contigs"].append(contig_ID)
                    if COG[0] != "-":
                        cluster_dict[cluster]["COGs"].extend(COG)


def output_writing(output, cluster_dict):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("[D] - Cell cycle control, cell division, chromosome partitioning\n"
                          "[M] - Cell wall/membrane/envelope biogenesis\n"
                          "[N] - Cell motility\n"
                          "[O] - Post-translational modification, protein turnover, and chaperones\n"
                          "[T] - Signal transduction mechanisms\n"
                          "[U] - Intracellular trafficking, secretion, and vesicular transport\n"
                          "[V] - Defense mechanisms\n"
                          "[W] - Extracellular structures\n"
                          "[Y] - Nuclear structure\n"
                          "[Z] - Cytoskeleton\n"
                          "[A] - RNA processing and modification\n"
                          "[B] - Chromatin structure and dynamics\n"
                          "[J] - Translation, ribosomal structure and biogenesis\n"
                          "[K] - Transcription\n"
                          "[L] - Replication, recombination and repair\n"
                          "[C] - Energy production and conversion\n"
                          "[E] - Amino acid transport and metabolism\n"
                          "[F] - Nucleotide transport and metabolism\n"
                          "[G] - Carbohydrate transport and metabolism\n"
                          "[H] - Coenzyme transport and metabolism\n"
                          "[I] - Lipid transport and metabolism\n"
                          "[P] - Inorganic ion transport and metabolism\n"
                          "[Q] - Secondary metabolites biosynthesis, transport, and catabolism\n"
                          "[R] - General function prediction only\n"
                          "[S] - Function unknown\n"
                          "Clusters\tD\tM\tN\tO\tT\tU\tV\tW\tY\tZ\tA\tB\tJ\tK\tL\tC\tE\tF\tG\tH\tI\tP\tQ\tR\tS\n")
        COGs = ["D", "M", "N", "O", "T", "U", "V", "W", "Y", "Z", "A", "B", "J", "K", "L", "C", "E", "F", "G", "H",
                "I", "P", "Q", "R", "S"]
        for cluster, values in cluster_dict.items():
            counts = ["{count}".format(count=values["COGs"].count(COG)) for COG in COGs]
            output_file.write("{ID}\t{counts}\n".format(ID=cluster, counts="\t".join(counts)))


if __name__ == "__main__":
    cluster_dict = {}
    print("***** Parsing input table *****")
    table_parsing(args.tab, cluster_dict)
    print("***** Output file creating *****")
    output_writing(args.output, cluster_dict)