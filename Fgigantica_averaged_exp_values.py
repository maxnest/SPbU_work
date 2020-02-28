try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--exp_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(gene_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split(";")
        #print(description)
        gene_id, metacercaria, egg, redia, cercaria, seventy_days_juvenile, miracidia, \
        marita, forty_two_days_juvenile = description[0], \
                np.mean([float(description[1]), float(description[2]), float(description[3])]), \
                np.mean([float(description[4]), float(description[5]), float(description[6])]), \
                np.mean([float(description[7]), float(description[8]), float(description[9])]), \
                np.mean([float(description[10]), float(description[11]), float(description[12])]), \
                np.mean([float(description[13]), float(description[14]), float(description[15])]), \
                np.mean([float(description[16]), float(description[17]), float(description[18])]), \
                np.mean([float(description[19]), float(description[20]), float(description[21])]), \
                np.mean([float(description[22]), float(description[23]), float(description[24])])
        gene_dict[gene_id] = {"metacercaria": metacercaria, "egg": egg, "redia": redia, "cercaria": cercaria,
                              "seventy_days": seventy_days_juvenile, "miracidia": miracidia, "marita": marita,
                              "forty_two_days": forty_two_days_juvenile}


def output_writing(gene_dict, out):
    with open("{out}.averaged_exp_values.tsv".format(out=out), 'a') as output:
        output.write("GeneID\tEgg\tMiracidia\tRedia\tCercaria\tMetacercaria\tForty_two_days_juvenile\t"
                     "Seventy_days_juvenile\tMarita\n")
        for gene, stages in gene_dict.items():
            output.write("{gene}\t{egg}\t{miracidia}\t{redia}\t{cercaria}\t{metacercaria}\t{forty_two}\t{seventy}\t"
                         "{marita}\n".format(gene=gene, egg=stages["egg"], miracidia=stages["miracidia"],
                                             redia=stages["redia"], cercaria=stages["cercaria"],
                                             metacercaria=stages["metacercaria"], forty_two=stages["forty_two_days"],
                                             seventy=stages["seventy_days"], marita=stages["marita"]))


if __name__ == "__main__":
    gene_dict = {}
    table_parsing(gene_dict, args.exp_table)
    output_writing(gene_dict, args.out)