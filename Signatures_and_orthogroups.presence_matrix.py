try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True,
                    help="The prepared table with orthogroups (OrthoFinder output). "
                         "The gene IDs should be in table "
                         "and species tag should be identical to tags in files with signatures "
                         "(for instance, 'Psim' for P.simillimum)")
parser.add_argument('--sign', type=argparse.FileType('r'), required=True,
                    help="The prepared table with molecular signatures (gene IDs and stage tags). "
                         "The stage tags should include species tag used in orthogroup table!"
                         "(for instance, 'Psim_redia' for redia of P. simillimum)")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def ortho_parsing(ortho, ortho_dict):
    header = ortho.readline()
    for line in ortho:
        description = line.strip().split("\t")
        OG, species_values = description[0], description[1:]
        for species_value in species_values:
            for geneID in species_value.split(", "):
                ortho_dict[geneID] = OG


def signature_parsing(sign, ortho_dict, sign_dict):
    header = sign.readline()
    for line in sign:
        geneID, stage_tag = line.strip().split("\t")[0], line.strip().split("\t")[1]
        if stage_tag not in sign_dict:
            sign_dict[stage_tag] = []

        if geneID in ortho_dict:
            sign_dict[stage_tag].append(ortho_dict[geneID])


def output_writing(output, sign_dict):
    all_stages = [stage for stage in sign_dict]
    all_ortho = []
    for value in sign_dict.values():
        all_ortho.extend(value)

    with open("{output}.presence_matrix.tsv".format(output=output), 'a') as output_file:
        output_file.write("Orthogroup\t{stages}\n".format(stages="\t".join(all_stages)))
        for og in set(all_ortho):
            values = ["1" if og in sign_dict[stage] else "0" for stage in all_stages]
            output_file.write("{og}\t{values}\n".format(og=og, values="\t".join(values)))


if __name__ == "__main__":
    ortho_dict, sign_dict = {}, {}
    print("***** Tab with orthogroups parsing *****")
    ortho_parsing(args.ortho, ortho_dict)
    print("***** Tab with signatures parsing ***** ")
    signature_parsing(args.sign, ortho_dict, sign_dict)
    print("***** Output file creating *****")
    output_writing(args.output, sign_dict)