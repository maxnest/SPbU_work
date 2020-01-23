try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--psilo_fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta file with Psilotrema simillimum nucleotide sequences")
parser.add_argument('--sphaer_fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta file with Sphaeridiotrema pseudoglobulus nucleotide sequences")
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="OrthoFinder output file")
parser.add_argument('--psilo_red_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of P.simillimum rediae-specific sequences")
parser.add_argument('--psilo_cer_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of P.simillimum cercariae-specific sequences")
parser.add_argument('--psilo_mar_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of P.simillimum marita-specific sequences")
parser.add_argument('--sphaer_red_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of S.pseudoglobulus rediae-specific sequences")
parser.add_argument('--sphaer_cer_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of S.pseudoglobulus cercariae-specific sequences")
parser.add_argument('--sphaer_mar_set', type=argparse.FileType('r'), required=True,
                    help="File with IDs of S.pseudoglobulus marita-specific sequences")
parser.add_argument('--psilo_red_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum rediae-specific gene set GO-enrichment analysis")
parser.add_argument('--psilo_cer_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum cercariae-specific gene set GO-enrichment analysis")
parser.add_argument('--psilo_mar_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of P.simillimum marita-specific gene set GO-enrichment analysis")
parser.add_argument('--sphaer_red_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.pseudoglobulus rediae-specific gene set GO-enrichment analysis")
parser.add_argument('--sphaer_cer_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.pseudoglobulus cercariae-specific gene set GO-enrichment analysis")
parser.add_argument('--sphaer_mar_go', type=argparse.FileType('r'), required=True,
                    help="Table with results of S.pseudoglobulus marita-specific gene set GO-enrichment analysis")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ortho": []}


def orthogroups_parsing(psilo_contig_dict, sphaer_contig_dict, orthogroups):
    head = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(","):
                if protein.split(".p")[0] in psilo_contig_dict.keys():
                    psilo_contig_dict[protein.split(".p")[0]]["ortho"].append(group_ID)
                elif protein.split(".p")[0] in sphaer_contig_dict.keys():
                    sphaer_contig_dict[protein.split(".p")[0]]["ortho"].append(group_ID)


def gene_set_parsing(file, list):
    for line in file:
        contig = line.strip()
        list.append(contig)


def go_table_parsing(file, list):
    head = file.readline()
    for line in file:
        description = line.strip().split("\t")
        list.append(description[0])


def jaccard_similarity(one_list, other_list):
    """
    FROM: http://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
    The function measure the similarity between two sets of genes (Jaccard similarity index)
    @param one_list: list of genes
    @param other_list: list of genes
    """
    intersection_cardinality = len(set.intersection(*[set(one_list), set(other_list)]))
    union_cardinality = len(set.union(*[set(one_list), set(other_list)]))
    if union_cardinality != 0:
        return round(intersection_cardinality / float(union_cardinality), 2)
    else:
        return 0


def gene_sets_comparison(psilo_ortho, sphaer_ortho, psilo_set, sphaer_set, psilo_tag, sphaer_tag, set_results_dict):
    psilo_ortho_in_set, sphaer_ortho_in_set = \
        set([psilo_ortho[contig]["ortho"][0] for contig in psilo_set if len(psilo_ortho[contig]["ortho"]) != 0]), \
        set([sphaer_ortho[contig]["ortho"][0] for contig in sphaer_set if len(sphaer_ortho[contig]["ortho"]) != 0])
    # set_results_dict["{psilo_tag}_vs_{sphaer_tag}".format(psilo_tag=psilo_tag, sphaer_tag=sphaer_tag)] = \
    #     jaccard_similarity(psilo_ortho_in_set, sphaer_ortho_in_set)
    set_results_dict["{psilo_tag}_vs_{sphaer_tag}".format(psilo_tag=psilo_tag, sphaer_tag=sphaer_tag)] = {
        "Jaccard": jaccard_similarity(psilo_ortho_in_set, sphaer_ortho_in_set),
        "Intersection": set.intersection(*[set(psilo_ortho_in_set), set(sphaer_ortho_in_set)])}


def go_enrichment_results_comparison(psilo_go, sphaer_go, psilo_tag, sphaer_tag, go_results_dict):
    # go_results_dict["{psilo_tag}_vs_{sphaer_tag}".format(psilo_tag=psilo_tag, sphaer_tag=sphaer_tag)] = \
    #     jaccard_similarity(psilo_go, sphaer_go)
    go_results_dict["{psilo_tag}_vs_{sphaer_tag}".format(psilo_tag=psilo_tag, sphaer_tag=sphaer_tag)] = {
        "Jaccard": jaccard_similarity(psilo_go, sphaer_go),
        "Intersection": set.intersection(*[set(psilo_go), set(sphaer_go)])}


def output_writing(output, result_dict, result_tag):
    with open("{output}.{result_tag}.tsv".format(output=output, result_tag=result_tag), 'a') as output_file:
        output_file.write("P.simillimum stages\S.pseudoglobulus stages\tRediae (Jaccard (Len of Intersection))\t"
                          "Cercariae (Jaccard (Len of Intersection))\tMarita (Jaccard (Len of Intersection))\n")
        output_file.write("Rediae\t{J_red_vs_red} ({I_red_vs_red})\t{J_red_vs_cer} ({I_red_vs_cer})\t"
                          "{J_red_vs_mar} ({I_red_vs_mar})\n"
                          "Cercariae\t{J_cer_vs_red} ({I_cer_vs_red})\t{J_cer_vs_cer} ({I_cer_vs_cer})\t"
                          "{J_cer_vs_mar} ({I_cer_vs_mar})\n"
                          "Marita\t{J_mar_vs_red} ({I_mar_vs_red})\t{J_mar_vs_cer} ({I_mar_vs_cer})\t"
                          "{J_mar_vs_mar} ({I_mar_vs_mar})".format(
                           J_red_vs_red=result_dict["psilo_red_vs_sphaer_red"]["Jaccard"],
                           I_red_vs_red=len(result_dict["psilo_red_vs_sphaer_red"]["Intersection"]),
                           J_red_vs_cer=result_dict["psilo_red_vs_sphaer_cer"]["Jaccard"],
                           I_red_vs_cer=len(result_dict["psilo_red_vs_sphaer_cer"]["Intersection"]),
                           J_red_vs_mar=result_dict["psilo_red_vs_sphaer_mar"]["Jaccard"],
                           I_red_vs_mar=len(result_dict["psilo_red_vs_sphaer_mar"]["Intersection"]),
                           J_cer_vs_red=result_dict["psilo_cer_vs_sphaer_red"]["Jaccard"],
                           I_cer_vs_red=len(result_dict["psilo_cer_vs_sphaer_red"]["Intersection"]),
                           J_cer_vs_cer=result_dict["psilo_cer_vs_sphaer_cer"]["Jaccard"],
                           I_cer_vs_cer=len(result_dict["psilo_cer_vs_sphaer_cer"]["Intersection"]),
                           J_cer_vs_mar=result_dict["psilo_cer_vs_sphaer_mar"]["Jaccard"],
                           I_cer_vs_mar=len(result_dict["psilo_cer_vs_sphaer_mar"]["Intersection"]),
                           J_mar_vs_red=result_dict["psilo_mar_vs_sphaer_red"]["Jaccard"],
                           I_mar_vs_red=len(result_dict["psilo_mar_vs_sphaer_red"]["Intersection"]),
                           J_mar_vs_cer=result_dict["psilo_mar_vs_sphaer_cer"]["Jaccard"],
                           I_mar_vs_cer=len(result_dict["psilo_mar_vs_sphaer_cer"]["Intersection"]),
                           J_mar_vs_mar=result_dict["psilo_mar_vs_sphaer_mar"]["Jaccard"],
                           I_mar_vs_mar=len(result_dict["psilo_mar_vs_sphaer_mar"]["Intersection"])))


if __name__ == "__main__":
    psilo_contig_dict, sphaer_contig_dict = {}, {}
    psilo_red_set, psilo_cer_set, psilo_mar_set, sphaer_red_set, sphaer_cer_set, sphaer_mar_set = [], [], [], [], [], []
    psilo_red_go, psilo_cer_go, psilo_mar_go, sphaer_red_go, sphaer_cer_go, sphaer_mar_go = [], [], [], [], [], []
    set_comparison_results, go_comparison_results = {}, {}
    print("***** Input files parsing *****")
    nucl_parsing(psilo_contig_dict, args.psilo_fasta)
    nucl_parsing(sphaer_contig_dict, args.sphaer_fasta)
    orthogroups_parsing(psilo_contig_dict, sphaer_contig_dict, args.ortho)
    gene_set_parsing(args.psilo_red_set, psilo_red_set)
    gene_set_parsing(args.psilo_cer_set, psilo_cer_set)
    gene_set_parsing(args.psilo_mar_set, psilo_mar_set)
    gene_set_parsing(args.sphaer_red_set, sphaer_red_set)
    gene_set_parsing(args.sphaer_cer_set, sphaer_cer_set)
    gene_set_parsing(args.sphaer_mar_set, sphaer_mar_set)
    go_table_parsing(args.psilo_red_go, psilo_red_go)
    go_table_parsing(args.psilo_cer_go, psilo_cer_go)
    go_table_parsing(args.psilo_mar_go, psilo_mar_go)
    go_table_parsing(args.sphaer_red_go, sphaer_red_go)
    go_table_parsing(args.sphaer_cer_go, sphaer_cer_go)
    go_table_parsing(args.sphaer_mar_go, sphaer_mar_go)
    print("***** Comparisons performing *****")
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_red_set, sphaer_red_set,
                         "psilo_red", "sphaer_red", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_red_set, sphaer_cer_set,
                         "psilo_red", "sphaer_cer", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_red_set, sphaer_mar_set,
                         "psilo_red", "sphaer_mar", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_cer_set, sphaer_red_set,
                         "psilo_cer", "sphaer_red", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_cer_set, sphaer_cer_set,
                         "psilo_cer", "sphaer_cer", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_cer_set, sphaer_mar_set,
                         "psilo_cer", "sphaer_mar", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_mar_set, sphaer_red_set,
                         "psilo_mar", "sphaer_red", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_mar_set, sphaer_cer_set,
                         "psilo_mar", "sphaer_cer", set_comparison_results)
    gene_sets_comparison(psilo_contig_dict, sphaer_contig_dict, psilo_mar_set, sphaer_mar_set,
                         "psilo_mar", "sphaer_mar", set_comparison_results)
    go_enrichment_results_comparison(psilo_red_go, sphaer_red_go, "psilo_red", "sphaer_red", go_comparison_results)
    go_enrichment_results_comparison(psilo_red_go, sphaer_cer_go, "psilo_red", "sphaer_cer", go_comparison_results)
    go_enrichment_results_comparison(psilo_red_go, sphaer_mar_go, "psilo_red", "sphaer_mar", go_comparison_results)
    go_enrichment_results_comparison(psilo_cer_go, sphaer_red_go, "psilo_cer", "sphaer_red", go_comparison_results)
    go_enrichment_results_comparison(psilo_cer_go, sphaer_cer_go, "psilo_cer", "sphaer_cer", go_comparison_results)
    go_enrichment_results_comparison(psilo_cer_go, sphaer_mar_go, "psilo_cer", "sphaer_mar", go_comparison_results)
    go_enrichment_results_comparison(psilo_mar_go, sphaer_red_go, "psilo_mar", "sphaer_red", go_comparison_results)
    go_enrichment_results_comparison(psilo_mar_go, sphaer_cer_go, "psilo_mar", "sphaer_cer", go_comparison_results)
    go_enrichment_results_comparison(psilo_mar_go, sphaer_mar_go, "psilo_mar", "sphaer_mar", go_comparison_results)
    print("***** Output creating *****")
    output_writing(args.output, set_comparison_results, "Jaccard_orthogroups_in_stage-specific_gene_sets")
    output_writing(args.output, go_comparison_results, "Jaccard_enriched_GO-terms_in_stage-specific_gene_sets")
