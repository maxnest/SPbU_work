try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--psilo_summary_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--sphaer_summary_table', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, red, cer, mar):
    for line in table:
        if not line.startswith("#") and not line.startswith("Contig_ID"):
            description = line.strip().split("\t")
            orthogroup, red_specific, cer_specific, mar_specific = description[1], description[-4], description[-3], \
                                                                   description[-2]
            if orthogroup != "-":
                if red_specific == "*":
                    red.append(orthogroup)

                if cer_specific == "*":
                    cer.append(orthogroup)

                if mar_specific == "*":
                    mar.append(orthogroup)


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


def orthosets_comparison(first_set, second_set, first_tag, second_tag, Jaccard_dict):
    intersection = set.intersection(*[set(first_set), set(second_set)])
    Jaccard_dict["{first_tag}_vs_{second_tag}".format(first_tag=first_tag, second_tag=second_tag)] = {
        "Jaccard": jaccard_similarity(first_set, second_set),
        "Intersection": intersection,
        "{first_tag}_specific".format(first_tag=first_tag): [term for term in first_set if term not in intersection],
        "{second_tag}_specific".format(second_tag=second_tag): [term for term in second_set if term not in intersection]}


def output_files_creating(Jaccard_dict, first_sp_tag, second_sp_tag, output):
    with open("{output}.tsv".format(output=output), 'a') as Jaccard_output:
        Jaccard_output.write("{first_tag}\{second_tag}\tRediae\tCercariae\tMarita\n".format(first_tag=first_sp_tag,
                                                                                            second_tag=second_sp_tag))
        Jaccard_output.write("Rediae\t{red_vs_red}\t{red_vs_cer}\t{red_vs_mar}\n".format(
                                red_vs_red=Jaccard_dict["{first}_rediae_vs_{second}_rediae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                red_vs_cer=Jaccard_dict["{first}_rediae_vs_{second}_cercariae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                red_vs_mar=Jaccard_dict["{first}_rediae_vs_{second}_marita".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"]))
        Jaccard_output.write("Cercariae\t{cer_vs_red}\t{cer_vs_cer}\t{cer_vs_mar}\n".format(
                                cer_vs_red=Jaccard_dict["{first}_cercariae_vs_{second}_rediae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                cer_vs_cer=Jaccard_dict["{first}_cercariae_vs_{second}_cercariae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                cer_vs_mar=Jaccard_dict["{first}_cercariae_vs_{second}_marita".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"]))
        Jaccard_output.write("Marita\t{mar_vs_red}\t{mar_vs_cer}\t{mar_vs_mar}\n".format(
                                mar_vs_red=Jaccard_dict["{first}_marita_vs_{second}_rediae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                mar_vs_cer=Jaccard_dict["{first}_marita_vs_{second}_cercariae".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"],
                                mar_vs_mar=Jaccard_dict["{first}_marita_vs_{second}_marita".format(
                                    first=first_sp_tag, second=second_sp_tag)]["Jaccard"]))

    for pair, values in Jaccard_dict.items():
        with open("{pair}.orthogroups_intersection.tsv".format(pair=pair), 'a') as pair_output:
            for orthogroup in values["Intersection"]:
                pair_output.write("{ortho}\n".format(ortho=orthogroup))


if __name__ == "__main__":
    Jaccard_dict, psilo_red, psilo_cer, psilo_mar, sphaer_red, sphaer_cer, sphaer_mar = {}, [], [], [], [], [], []
    table_parsing(args.psilo_summary_table, psilo_red, psilo_cer, psilo_mar)
    table_parsing(args.sphaer_summary_table, sphaer_red, sphaer_cer, sphaer_mar)
    orthosets_comparison(psilo_red, sphaer_red, "Psimillimum_rediae", "Spseudoglobulus_rediae", Jaccard_dict)
    orthosets_comparison(psilo_red, sphaer_cer, "Psimillimum_rediae", "Spseudoglobulus_cercariae", Jaccard_dict)
    orthosets_comparison(psilo_red, sphaer_mar, "Psimillimum_rediae", "Spseudoglobulus_marita", Jaccard_dict)
    orthosets_comparison(psilo_cer, sphaer_red, "Psimillimum_cercariae", "Spseudoglobulus_rediae", Jaccard_dict)
    orthosets_comparison(psilo_cer, sphaer_cer, "Psimillimum_cercariae", "Spseudoglobulus_cercariae", Jaccard_dict)
    orthosets_comparison(psilo_cer, sphaer_mar, "Psimillimum_cercariae", "Spseudoglobulus_marita", Jaccard_dict)
    orthosets_comparison(psilo_mar, sphaer_red, "Psimillimum_marita", "Spseudoglobulus_rediae", Jaccard_dict)
    orthosets_comparison(psilo_mar, sphaer_cer, "Psimillimum_marita", "Spseudoglobulus_cercariae", Jaccard_dict)
    orthosets_comparison(psilo_mar, sphaer_mar, "Psimillimum_marita", "Spseudoglobulus_marita", Jaccard_dict)
    output_files_creating(Jaccard_dict, "Psimillimum", "Spseudoglobulus", args.out)