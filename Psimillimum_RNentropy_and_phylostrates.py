try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--psim_red', type=argparse.FileType('r'), required=True,
                    help="Table with sequences over-expressed in P.simillimum redia stage (RNentropy results)")
parser.add_argument('--psim_cer', type=argparse.FileType('r'), required=True,
                    help="Table with sequences over-expressed in P.simillimum cercaria stage (RNentropy results)")
parser.add_argument('--psim_mar', type=argparse.FileType('r'), required=True,
                    help="Table with sequences over-expressed in P.simillimum adult worm stage (RNentropy results)")
parser.add_argument('--gene_map', type=argparse.FileType('r'), required=True,
                    help="Gene map: GeneID, TranscriptID, ProteinID")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum phylostratigraphy analysis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def RNentropy_parsing(table, list):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene_id, values = description[0][1:-1], description[1:]
        list.append(gene_id)


def gene_map_parsing(gene_map, gene_dict):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_id, transcript_id, protein_id = description[0], description[1], description[2]
        gene_dict[gene_id] = protein_id


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        phylostratr_dict[protein_ID] = mrca_name


def append_gene(list, merged_dict, phylostratr_dict, gene_dict, tag):
    for gene_id in list:
        merged_dict[phylostratr_dict[gene_dict[gene_id]]][tag].append(gene_id)


def results_merging(phylostratr_dict, gene_dict, redia_list, cercaria_list, marita_list, merged_dict):
    for protein_ID, mrca_name in phylostratr_dict.items():
        if mrca_name not in merged_dict.keys():
            merged_dict[mrca_name] = {"Redia": [], "Cercaria": [], "Adult_worm": []}

    append_gene(redia_list, merged_dict, phylostratr_dict, gene_dict, "Redia")
    append_gene(cercaria_list, merged_dict, phylostratr_dict, gene_dict, "Cercaria")
    append_gene(marita_list, merged_dict, phylostratr_dict, gene_dict, "Adult_worm")

    for key, values in merged_dict.items():
        for stage, genes in values.items():
            if len(genes) == 0:
                genes.append("-")


def output_writing(output, merged_dict, redia_list, cercaria_list, marita_list):
    with open("{output}.overexpressed_in_phylostrates.tsv".format(output=output), "a") as output_file:
        output_file.write("Phylostrates\tRedia:gene_count\tRedia:percent\tCercaria:gene_count\tCercaria:percent\t"
                          "Adult_worm:gene_count\tAdult_worm:percent\n")
        for mrca_name, values in merged_dict.items():
            output_file.write("{mrca_name}\t{red_count}\t{red_percent}\t{cer_count}\t{cer_percent}\t"
                              "{adult_count}\t{adult_percent}\n".format(
                                mrca_name=mrca_name,
                                red_count=len(values["Redia"]),
                                red_percent=round((len(values["Redia"])/len(redia_list))*100, 2),
                                cer_count=len(values["Cercaria"]),
                                cer_percent=round((len(values["Cercaria"])/len(cercaria_list))*100, 2),
                                adult_count=len(values["Adult_worm"]),
                                adult_percent=round((len(values["Adult_worm"])/len(marita_list))*100, 2)))


if __name__ == "__main__":
    redia_list, cercaria_list, marita_list = [], [], []
    gene_dict, phylostratr_dict, merged_dict = {}, {}, {}
    print("***** Input tables parsing *****")
    RNentropy_parsing(args.psim_red, redia_list)
    RNentropy_parsing(args.psim_cer, cercaria_list)
    RNentropy_parsing(args.psim_mar, marita_list)
    gene_map_parsing(args.gene_map, gene_dict)
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    print("***** Results merging *****")
    results_merging(phylostratr_dict, gene_dict, redia_list, cercaria_list, marita_list, merged_dict)
    print("***** Output file writing *****")
    output_writing(args.output, merged_dict, redia_list, cercaria_list, marita_list)