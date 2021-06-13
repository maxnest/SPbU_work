try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--gene_map', type=argparse.FileType('r'), required=True,
                    help="Table with 3 columns: gene_ID, transcript_ID, protein_ID")
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True)
parser.add_argument('--ortho_anno', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nt', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nr', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--expression', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_overexp', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_overexp', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_overexp', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_overexp', type=argparse.FileType('r'), required=True)
parser.add_argument('--classic_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--nonclassic_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def gene_map_parsing(gene_dict, gene_map):
    header = gene_map.readline()
    for line in gene_map:
        description = line.strip().split("\t")
        gene_ID, transcript_ID, protein_ID = description[0], description[1], description[2]
        gene_dict[gene_ID] = {"transcript": transcript_ID, "protein": protein_ID,
                              "domains": [],
                              "Anno_swiss": {"hit": [], "identity": []},
                              "Anno_nt": {"hit": [], "identity": []},
                              "Anno_nr": {"hit": [], "identity": []},
                              "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                              "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                              "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                              "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                              "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                              "Externa_TPMs": 0, "Growing_TPMs": 0, "Middle_TPMs": 0,
                              "Terminal_TPMs": 0, "Over-expression": [],
                              "Ortho": {"id": [],
                                        "anno": {"GO": {"ids": [], "desc": []},
                                                 "PfamA": {"ids": [], "desc": []},
                                                 "IPR": {"ids": [], "desc": []}}},
                              "ExSec": [],
                              "Phylostratum": []}


def orthogroups_parsing(gene_dict, prot_2_gene_dict, orthogroups):
    header = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(", "):
                if protein in prot_2_gene_dict.keys():
                    gene_dict[prot_2_gene_dict[protein]]["Ortho"]["id"].append(group_ID)

    for gene, values in gene_dict.items():
        if len(values["Ortho"]["id"]) == 0:
            values["Ortho"]["id"].append("-")


def orthogroups_annotation(gene_dict, ortho_anno):
    header = ortho_anno.readline()
    for line in ortho_anno:
        description = line.strip().split("\t")
        ortho_id, go_ids, go_desc, pfam_ids, pfam_desc, ipr_ids, ipr_desc = \
            description[0][1:-1], description[3][1:-1], description[4][1:-1], description[5][1:-1], \
            description[6][1:-1], description[7][1:-1], description[8][1:-1]
        for gene, values in gene_dict.items():
            if values["Ortho"]["id"][0] == ortho_id:
                values["Ortho"]["anno"]["GO"]["ids"].append(go_ids)
                values["Ortho"]["anno"]["GO"]["desc"].append(go_desc)
                values["Ortho"]["anno"]["PfamA"]["ids"].append(pfam_ids)
                values["Ortho"]["anno"]["PfamA"]["desc"].append(pfam_desc)
                values["Ortho"]["anno"]["IPR"]["ids"].append(ipr_ids)
                values["Ortho"]["anno"]["IPR"]["desc"].append(ipr_desc)
                # print(values["Ortho"]["anno"])

    for gene, values in gene_dict.items():
        if values["Ortho"]["id"][0] == "-":
            values["Ortho"]["anno"]["GO"]["ids"].append("-")
            values["Ortho"]["anno"]["GO"]["desc"].append("-")
            values["Ortho"]["anno"]["PfamA"]["ids"].append("-")
            values["Ortho"]["anno"]["PfamA"]["desc"].append("-")
            values["Ortho"]["anno"]["IPR"]["ids"].append("-")
            values["Ortho"]["anno"]["IPR"]["desc"].append("-")


def phylostratr_parsing(gene_dict, prot_2_gene_dict, phylostratr):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        protein_ID, mrca, ps, mrca_name = description[0][1:-1], description[1], description[2], description[3][1:-1]
        if protein_ID in prot_2_gene_dict.keys():
            gene_dict[prot_2_gene_dict[protein_ID]]["Phylostratum"].append("{ps}:{mrca_name}".format(
                ps=ps, mrca_name=mrca_name))

    for gene, values in gene_dict.items():
        if len(values["Phylostratum"]) == 0:
            values["Phylostratum"].append("-")


def eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, eggNOG):
    keys = ["EggNOG:Preferred_name", "EggNOG:GO_terms", "EggNOG:EC_number", "EggNOG:KEGG_KO", "EggNOG:KEGG_Pathway",
            "EggNOG:KEGG_Module", "EggNOG:KEGG_Reaction", "EggNOG:rclass", "EggNOG:BRITE", "EggNOG:KEGG_TC",
            "EggNOG:CAZy", "EggNOG:BiGG_Reaction", "EggNOG:OG", "EggNOG:COG_cat", "EggNOG:Description"]
    index_in_annotation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20]

    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            protein, annotation = description[0], description[1:]
            for el in index_in_annotation:
                if len(annotation) >= el + 1 and len(annotation[el]) != 0:
                    if protein in prot_2_gene_dict.keys():
                        gene_dict[prot_2_gene_dict[protein]][keys[index_in_annotation.index(el)]].append(annotation[el])

    for gene, values in gene_dict.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append('-')


def add_match(description, ID, hit_name, gene_dict, seq_dict, key_tag):
    if hit_name == "no hits":
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["hit"].append("-")
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["identity"].append("-")
    else:
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
        gene_dict[seq_dict[ID]]["Anno_{key}".format(key=key_tag)]["identity"].append(description[9])


def BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, BLAST, key_tag):
    header = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        ID, hit_name = description[0], description[3]
        if ID in trans_2_gene_dict.keys():
            add_match(description, ID, hit_name, gene_dict, trans_2_gene_dict, key_tag)
        elif ID in prot_2_gene_dict.keys():
            add_match(description, ID, hit_name, gene_dict, prot_2_gene_dict, key_tag)


def domains_parsing(gene_dict, prot_2_gene_dict, domains):
    header = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        protein, domain_name = description[0].strip(), description[1]
        if protein in prot_2_gene_dict.keys():
            gene_dict[prot_2_gene_dict[protein]]["domains"].append(domain_name)

    for gene in gene_dict.keys():
        if len(gene_dict[gene]["domains"]) == 0:
            gene_dict[gene]["domains"].append("-")


def RNentropy_parsing(gene_dict, table, tag):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene = description[0][1:-1]
        gene_dict[gene]["Over-expression"].append(tag)


def expression_parsing(gene_dict, table):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        gene, externa, growing, middle, terminal, whole_body = \
            description[0], float(description[1]), float(description[2]), float(description[3]), \
            float(description[4]), float(description[5])
        gene_dict[gene]["Externa_TPMs"] += externa
        gene_dict[gene]["Growing_TPMs"] += growing
        gene_dict[gene]["Middle_TPMs"] += middle
        gene_dict[gene]["Terminal_TPMs"] += terminal


def exsec_parsing(gene_dict, exsec_seqs, prot_2_gene_dict, tag):
    for fasta in exsec_seqs:
        description = fasta.id.strip().split(" ")
        gene_dict[prot_2_gene_dict[description[0]]]["ExSec"].append(tag)


def write_output_files(gene_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("Gene_ID\tTranscript_ID\tProtein_ID\tPhylostrata\tOrthogroup\tOrthogroup:GOterms_IDs\t"
                          "Orthogroup:GOterms_description\tOrthogroup:PfamA_domains\tOrthogroup:PfamA_description\t"
                          "Orthogroup:IPR_domains\tOrthogroup:IPR_description\t"
                          "NCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\t"
                          "Domains_arch\tEggNOG:Preferred_name\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tExterna_averaged_TPMs\t"
                          "Growing_trunk_part_averaged_TPMs\tMain_trunk_part_averaged_TPMs\t"
                          "Thoracic_part_of_interna_averaged_TPMs\t"
                          "Over-expression\tSecretory\Excretory_sequences\n")
        for gene, values in gene_dict.items():
            output_file.write("{gene}\t{trans}\t{prot}\t{phylo}\t{ortho}\t{ortho_go}\t{ortho_go_desc}\t"
                              "{ortho_pfam}\t{ortho_pfam_desc}\t{ortho_ipr}\t{ortho_ipr_desc}\t"
                              "{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t"
                              "{swiss}\t{swiss_identity}\t{domains}\t{name}\t{go}\t"
                              "{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t"
                              "{og}\t{cog}\t{description}\t{externa}\t{growing}\t{middle}\t{terminal}\t"
                              "{overexp}\t{exsec}\n".format(
                                gene=gene, trans=values["transcript"], prot=values["protein"],
                                phylo=values["Phylostratum"][0],
                                ortho=values["Ortho"]["id"][0],
                                ortho_go=values["Ortho"]["anno"]["GO"]["ids"][0],
                                ortho_go_desc=values["Ortho"]["anno"]["GO"]["desc"][0],
                                ortho_pfam=values["Ortho"]["anno"]["PfamA"]["ids"][0],
                                ortho_pfam_desc=values["Ortho"]["anno"]["PfamA"]["desc"][0],
                                ortho_ipr=values["Ortho"]["anno"]["IPR"]["ids"][0],
                                ortho_ipr_desc=values["Ortho"]["anno"]["IPR"]["desc"][0],
                                nt=values["Anno_nt"]["hit"][0], nt_identity=values["Anno_nt"]["identity"][0],
                                nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                                swiss=values["Anno_swiss"]["hit"][0],
                                swiss_identity=values["Anno_swiss"]["identity"][0],
                                domains="|".join(values["domains"]),
                                name=values["EggNOG:Preferred_name"][0], go=values["EggNOG:GO_terms"][0],
                                ec=values["EggNOG:EC_number"][0], ko=values["EggNOG:KEGG_KO"][0],
                                pathway=values["EggNOG:KEGG_Pathway"][0], module=values["EggNOG:KEGG_Module"][0],
                                reaction=values["EggNOG:KEGG_Reaction"][0], rclass=values["EggNOG:rclass"][0],
                                brite=values["EggNOG:BRITE"][0], tc=values["EggNOG:KEGG_TC"][0],
                                cazy=values["EggNOG:CAZy"][0], bigg=values["EggNOG:BiGG_Reaction"][0],
                                og=values["EggNOG:OG"][0], cog=values["EggNOG:COG_cat"][0],
                                description=values["EggNOG:Description"][0], externa=values["Externa_TPMs"],
                                growing=values["Growing_TPMs"], middle=values["Middle_TPMs"],
                                terminal=values["Terminal_TPMs"],
                                overexp="|".join(values["Over-expression"]), exsec=values["ExSec"][0]))


if __name__ == "__main__":
    gene_dict, cluster_dict = {}, {}
    gene_map_parsing(gene_dict, args.gene_map)
    prot_2_gene_dict = {gene_dict[gene]["protein"]: gene for gene in gene_dict.keys()}
    trans_2_gene_dict = {gene_dict[gene]["transcript"]: gene for gene in gene_dict.keys()}
    print("***** Orthogroups parsing *****")
    orthogroups_parsing(gene_dict, prot_2_gene_dict, args.ortho)
    print("***** Orthogroups annotation parsing *****")
    orthogroups_annotation(gene_dict, args.ortho_anno)
    print("***** Phylostratigraphy results parsing *****")
    phylostratr_parsing(gene_dict, prot_2_gene_dict, args.phylostratr)
    print("***** eggNOG-mapper output parsing *****")
    eggNOG_mapper_parsing(gene_dict, prot_2_gene_dict, args.eggnog)
    print("***** BLAST results parsing *****")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_swiss, "swiss")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nt, "nt")
    BLAST_annotation(gene_dict, prot_2_gene_dict, trans_2_gene_dict, args.blast_nr, "nr")
    print("***** Domain architecture parsing *****")
    domains_parsing(gene_dict, prot_2_gene_dict, args.domains)
    print("***** Files with IDs of sequences with preferential expression patterns parsing *****")
    RNentropy_parsing(gene_dict, args.externa_overexp, "Externa")
    RNentropy_parsing(gene_dict, args.growing_overexp, "Growing_trunk")
    RNentropy_parsing(gene_dict, args.middle_overexp, "Main_trunk_part")
    RNentropy_parsing(gene_dict, args.terminal_overexp, "Thoracic_part_of_interna")
    for gene, values in gene_dict.items():
        if len(values["Over-expression"]) == 0:
            values["Over-expression"].append('-')
    print("***** Table with averaged expression values parsing *****")
    expression_parsing(gene_dict, args.expression)
    print("***** Parsing of fasta-files with excretory/secretory sequences *****")
    classic_seqs = SeqIO.parse(args.classic_exsec, "fasta")
    nonclassic_seqs = SeqIO.parse(args.nonclassic_exsec, "fasta")
    exsec_parsing(gene_dict, classic_seqs, prot_2_gene_dict, "Potential_classical")
    exsec_parsing(gene_dict, nonclassic_seqs, prot_2_gene_dict, "Potential_nonclassical")
    for gene, values in gene_dict.items():
        if len(values["ExSec"]) == 0:
            values["ExSec"].append('-')
    print("***** Output file creating *****")
    write_output_files(gene_dict, args.output)


