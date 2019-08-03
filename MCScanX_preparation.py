try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--orthogroups', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp1gff3', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp2gff3', type=argparse.FileType('r'), required=True)
parser.add_argument('--sp1tag', type=str, required=True)
parser.add_argument('--sp2tag', type=str, required=True)
args = parser.parse_args()


def gff_parsing(gff3, sp_tag, sp_dict):
    for line in gff3:
        if not line.startswith("#"):
            description = line.strip().split("\t") # tab-separated?
            if description[-1].startswith("ID=gene:"):
                scf, start, end, name = description[0], description[3], description[4], \
                                    description[-1].split(";")[0].split(":")[1]
                sp_dict[name] = {"scf": "{sp}|{scf}".format(sp=sp_tag, scf=scf), "start": start, "end": end}


def orthogroups_separation(orthogroups, orthodict, sp1_dict, sp2_dict, sp1_tag, sp2_tag):
    for line in orthogroups:
        description = line.strip().split(":")
        group_ID, genes = description[0], description[1]
        orthodict[group_ID] = {"{sp1}".format(sp1=sp1_tag): [], "{sp2}".format(sp2=sp2_tag): []}
        for gene in genes.split(" "):
            if gene in sp1_dict.keys():
                orthodict[group_ID]["{sp1}".format(sp1=sp1_tag)].append(gene)
            elif gene in sp2_dict.keys():
                orthodict[group_ID]["{sp2}".format(sp2=sp2_tag)].append(gene)


def pairs_searching(orthodict, orthopairs, sp1_tag, sp2_tag):
    for group_ID, values in orthodict.items():
        if len(values["{sp1}".format(sp1=sp1_tag)]) != 0 and len(values["{sp2}".format(sp2=sp2_tag)]) != 0:
            combinations = list(itertools.product(
                values["{sp1}".format(sp1=sp1_tag)], values["{sp2}".format(sp2=sp2_tag)]))
            for combination in combinations:
                combination_list = list(combination)
                orthopairs["{sp1_gene}_{sp2_gene}".format(
                    sp1_gene=combination_list[0], sp2_gene=combination_list[1])] = {
                    "{sp1}".format(sp1=sp1_tag): combination_list[0], "{sp2}".format(sp2=sp2_tag): combination_list[1]}


def write_homology(orthopairs, sp1_tag, sp2_tag):
    with open("{sp1}_{sp2}.homology".format(sp1=sp1_tag, sp2=sp2_tag), 'a') as homology:
        for key, value in orthopairs.items():
            homology.write("{sp1_gene}\t{sp2_gene}\n".format(
                sp1_gene=value["{sp1}".format(sp1=sp1_tag)], sp2_gene=value["{sp2}".format(sp2=sp2_tag)]))


def write_gff(orthopairs, sp1_dict, sp2_dict, sp1_tag, sp2_tag):
    sp1_genes, sp2_genes = [], []
    with open("{sp1}_{sp2}.gff".format(sp1=sp1_tag, sp2=sp2_tag), 'a') as gff:
        for key, value in orthopairs.items():
            sp1_gene = value["{sp1}".format(sp1=sp1_tag)]
            if sp1_gene not in sp1_genes:
                sp1_genes.append(sp1_gene)
                gff.write("{sp1_scf}\t{sp1_gene}\t{start}\t{end}\n".format(
                        sp1_scf=sp1_dict[sp1_gene]["scf"], sp1_gene=sp1_gene,
                        start=sp1_dict[sp1_gene]["start"], end=sp1_dict[sp1_gene]["end"]))

        for key, value in orthopairs.items():
            sp2_gene = value["{sp2}".format(sp2=sp2_tag)]
            if sp2_gene not in sp2_genes:
                sp2_genes.append(sp2_gene)
                gff.write("{sp2_scf}\t{sp2_gene}\t{start}\t{end}\n".format(
                            sp2_scf=sp2_dict[sp2_gene]["scf"], sp2_gene=sp2_gene,
                            start=sp2_dict[sp2_gene]["start"], end=sp2_dict[sp2_gene]["end"]))


if __name__ == "__main__":
    sp1_gff, sp2_gff, orthodict, orthopairs = {}, {}, {}, {}
    print("*** GFF3 parsing ***")
    gff_parsing(args.sp1gff3, args.sp1tag, sp1_gff)
    gff_parsing(args.sp2gff3, args.sp2tag, sp2_gff)
    print("*** Orthogroups separation ***")
    orthogroups_separation(args.orthogroups, orthodict, sp1_gff, sp2_gff, args.sp1tag, args.sp2tag)
    print("*** Orthopairs creating ***")
    pairs_searching(orthodict, orthopairs, args.sp1tag, args.sp2tag)
    print("*** Output files creating ***")
    write_homology(orthopairs, args.sp1tag, args.sp2tag)
    write_gff(orthopairs, sp1_gff, sp2_gff, args.sp1tag, args.sp2tag)

