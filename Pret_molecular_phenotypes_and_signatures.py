try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--norm', type=argparse.FileType('r'), required=True,
                    help="Table with normalized (rlog) mean expression values")
parser.add_argument('--whole_vs_externa_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_whole_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_externa_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_vs_whole_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_vs_externa_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_vs_whole_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_vs_externa_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_whole_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def molecular_phenotype(norm, phenotypes_dict):
    header = norm.readline()
    for line in norm:
        description = line.strip().split("\t")
        contig_ID, externa, growing, middle, terminal, whole_body = description[0], float(description[1]), \
                                                                    float(description[2]), float(description[3]), \
                                                                    float(description[4]), float(description[5])
        if externa > 0:
            phenotypes_dict["externa"][contig_ID] = externa

        if growing > 0:
            phenotypes_dict["growing"][contig_ID] = growing

        if middle > 0:
            phenotypes_dict["middle"][contig_ID] = middle

        if terminal > 0:
            phenotypes_dict["terminal"][contig_ID] = terminal

        if whole_body > 0:
            phenotypes_dict["whole_body"][contig_ID] = whole_body


def sign_parsing(dict, pair):
    header = pair.readline()
    for line in pair:
        description = line.strip().split("\t")
        contig_ID, LFC = description[0], float(description[2])
        dict[contig_ID] = LFC


def lfc_appending(sample_dict, key, pair_dict, tag):
    if key in pair_dict.keys():
        sample_dict[key][tag].append(pair_dict[key])
    else:
        sample_dict[key][tag].append("-")


def contrast_results_merging(sample_dict, first_dict, second_dict, third_dict, fourth_dict,
                             first_tag, second_tag, third_tag, fourth_tag):
    keys = []
    keys.extend(first_dict.keys())
    keys.extend(second_dict.keys())
    keys.extend(third_dict.keys())
    keys.extend(fourth_dict.keys())
    for key in set(keys):
        sample_dict[key] = {first_tag: [], second_tag: [], third_tag: [], fourth_tag: []}
        lfc_appending(sample_dict, key, first_dict, first_tag)
        lfc_appending(sample_dict, key, second_dict, second_tag)
        lfc_appending(sample_dict, key, third_dict, third_tag)
        lfc_appending(sample_dict, key, fourth_dict, fourth_tag)


def classification(sample_dict, signature_dict, target, first_tag, second_tag, third_tag, fourth_tag):
    for contig, values in sample_dict.items():
        if values[first_tag][0] != "-" and values[second_tag][0] != "-" and values[third_tag][0] != "-" and values[fourth_tag][0] != "-":
            if values[first_tag][0] > 0 and values[second_tag][0] > 0 and values[third_tag][0] > 0 and values[fourth_tag][0] > 0:
                signature_dict["core"].append(contig)
        else:
            if values[first_tag][0] != "-" and values[first_tag][0] > 0:
                signature_dict[first_tag].append(contig)

            if values[second_tag][0] != "-" and values[second_tag][0] > 0:
                signature_dict[second_tag].append(contig)

            if values[third_tag][0] != "-" and values[third_tag][0] > 0:
                signature_dict[third_tag].append(contig)

            if values[fourth_tag][0] != "-" and values[fourth_tag][0] > 0:
                signature_dict[fourth_tag].append(contig)

    print("***** {target} molecular signatures: ******\n"
          "Core: {core_count}\n"
          "{first}: {first_count}\n"
          "{second}: {second_count}\n"
          "{third}: {third_count}\n"
          "{fourth}: {fourth_count}\n".format(
           target=target, core_count=len(signature_dict["core"]),
           first=first_tag, first_count=len(signature_dict[first_tag]),
           second=second_tag, second_count=len(signature_dict[second_tag]),
           third=third_tag, third_count=len(signature_dict[third_tag]),
           fourth=fourth_tag, fourth_count=len(signature_dict[fourth_tag])))


def phenotype_writing(out, dict):
    for key, values in dict.items():
        with open("{out}.molecular_phenotype_of_{key}.tsv".format(out=out, key=key), 'a') as output:
            output.write("Contig_ID\tExpression(DESeq2:rlog normalization)\n")
            for contig, expression in values.items():
                output.write("{contig}\t{exp}\n".format(contig=contig, exp=expression))


def signatures_writing(out, signature_dict, target, first_tag, second_tag, third_tag, fourth_tag):
    with open("{out}.molecular_signature_of_{target}.tsv".format(out=out, target=target), 'a') as output:
        output.write("Core\t{core}\n{first_tag}\t{first}\n{second_tag}\t{second}\n{third_tag}\t{third}\n"
                     "{fourth_tag}\t{fourth}\n".format(core="\t".join(signature_dict["core"]),
                      first_tag=first_tag, first="\t".join(signature_dict[first_tag]),
                      second_tag=second_tag, second="\t".join(signature_dict[second_tag]),
                      third_tag=third_tag, third="\t".join(signature_dict[third_tag]),
                      fourth_tag=fourth_tag, fourth="\t".join(signature_dict[fourth_tag])))


if __name__ == "__main__":
    phenotypes_dict = {"externa": {}, "growing": {}, "middle": {}, "terminal": {}, "whole_body": {}}
    whole_vs_externa_sign, whole_vs_growing_sign, whole_vs_middle_sign, whole_vs_terminal_sign, \
        growing_vs_whole_sign, growing_vs_externa_sign, growing_vs_middle_sign, growing_vs_terminal_sign, \
        middle_vs_whole_sign, middle_vs_externa_sign, middle_vs_growing_sign, middle_vs_terminal_sign, \
        terminal_vs_whole_sign, terminal_vs_externa_sign, terminal_vs_growing_sign, terminal_vs_middle_sign, \
        externa_vs_whole_sign, externa_vs_growing_sign, externa_vs_middle_sign, externa_vs_terminal_sign = \
        {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    whole_body_dict, growing_dict, middle_dict, terminal_dict, externa_dict = {}, {}, {}, {}, {}
    whole_body_signature, growing_signature, middle_signature, terminal_signature, externa_signature = \
        {"core": [], "externa":[], "growing": [], "middle": [], "terminal": []}, \
        {"core": [], "whole_body": [], "externa": [], "middle": [], "terminal": []}, \
        {"core": [], "whole_body": [], "externa": [], "growing": [], "terminal": []}, \
        {"core": [], "whole_body": [], "externa": [], "growing": [], "middle": []}, \
        {"core": [], "whole_body": [], "growing": [], "middle": [], "terminal": []}
    print("***** Input tables parsing *****")
    molecular_phenotype(args.norm, phenotypes_dict)
    sign_parsing(whole_vs_externa_sign, args.whole_vs_externa_sign)
    sign_parsing(whole_vs_growing_sign, args.whole_vs_growing_sign)
    sign_parsing(whole_vs_middle_sign, args.whole_vs_middle_sign)
    sign_parsing(whole_vs_terminal_sign, args.whole_vs_terminal_sign)
    sign_parsing(growing_vs_whole_sign, args.growing_vs_whole_sign)
    sign_parsing(growing_vs_externa_sign, args.growing_vs_externa_sign)
    sign_parsing(growing_vs_middle_sign, args.growing_vs_middle_sign)
    sign_parsing(growing_vs_terminal_sign, args.growing_vs_terminal_sign)
    sign_parsing(middle_vs_whole_sign, args.middle_vs_whole_sign)
    sign_parsing(middle_vs_externa_sign, args.middle_vs_externa_sign)
    sign_parsing(middle_vs_growing_sign, args.middle_vs_growing_sign)
    sign_parsing(middle_vs_terminal_sign, args.middle_vs_terminal_sign)
    sign_parsing(terminal_vs_whole_sign, args.terminal_vs_whole_sign)
    sign_parsing(terminal_vs_externa_sign, args.terminal_vs_externa_sign)
    sign_parsing(terminal_vs_growing_sign, args.terminal_vs_growing_sign)
    sign_parsing(terminal_vs_middle_sign, args.terminal_vs_middle_sign)
    sign_parsing(externa_vs_whole_sign, args.externa_vs_whole_sign)
    sign_parsing(externa_vs_growing_sign, args.externa_vs_growing_sign)
    sign_parsing(externa_vs_middle_sign, args.externa_vs_middle_sign)
    sign_parsing(externa_vs_terminal_sign, args.externa_vs_terminal_sign)
    print("***** Results merging *****")
    contrast_results_merging(whole_body_dict, whole_vs_externa_sign, whole_vs_growing_sign, whole_vs_middle_sign,
                             whole_vs_terminal_sign, "externa", "growing", "middle", "terminal")
    contrast_results_merging(growing_dict, growing_vs_whole_sign, growing_vs_externa_sign, growing_vs_middle_sign,
                             growing_vs_terminal_sign, "whole_body", "externa", "middle", "terminal")
    contrast_results_merging(middle_dict, middle_vs_whole_sign, middle_vs_externa_sign, middle_vs_growing_sign,
                             middle_vs_terminal_sign, "whole_body", "externa", "growing", "terminal")
    contrast_results_merging(terminal_dict, terminal_vs_whole_sign, terminal_vs_externa_sign, terminal_vs_growing_sign,
                             terminal_vs_middle_sign, "whole_body", "externa", "growing", "middle")
    contrast_results_merging(externa_dict, externa_vs_whole_sign, externa_vs_growing_sign, externa_vs_middle_sign,
                             externa_vs_terminal_sign, "whole_body", "growing", "middle", "terminal")
    print("***** Sequence classification *****")
    classification(whole_body_dict, whole_body_signature, "Whole_body", "externa", "growing", "middle", "terminal")
    classification(growing_dict, growing_signature, "Growing_stolon", "whole_body", "externa", "middle", "terminal")
    classification(middle_dict, middle_signature, "Middle_part", "whole_body", "externa", "growing", "terminal")
    classification(terminal_dict, terminal_signature, "Terminal_part", "whole_body", "externa", "growing", "middle")
    classification(externa_dict, externa_signature, "Externa", "whole_body", "growing", "middle", "terminal")
    print("***** Output files creating *****")
    phenotype_writing(args.out, phenotypes_dict)
    signatures_writing(args.out, whole_body_signature, "whole_body", "externa", "growing", "middle", "terminal")
    signatures_writing(args.out, growing_signature, "growing_stolon", "whole_body", "externa", "middle", "terminal")
    signatures_writing(args.out, middle_signature, "middle_part", "whole_body", "externa", "growing", "terminal")
    signatures_writing(args.out, terminal_signature, "terminal_part", "whole_body", "externa", "growing", "middle")
    signatures_writing(args.out, externa_signature, "externa", "whole_body", "growing", "middle", "terminal")

