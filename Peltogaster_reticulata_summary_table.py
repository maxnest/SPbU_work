try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--class_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--nonclass_exsec', type=argparse.FileType('r'), required=True)
parser.add_argument('--ortho', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nt', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nr', type=argparse.FileType('r'), required=True)
parser.add_argument('--neuropep', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_jongeneel', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_jongeneel', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_jongeneel', type=argparse.FileType('r'), required=True)
parser.add_argument('--terminal_jongeneel', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_jongeneel', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_externa_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--whole_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--growing_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--middle_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_growing_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_middle_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--externa_vs_terminal_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ID": [seq.id.split(" ")[0]], "nucl": seq.seq, "amino": "",
                               "domains": [],  "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [],
                               "EggNOG:EC_number": [], "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [],
                               "EggNOG:KEGG_Module": [], "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [],
                               "EggNOG:BRITE": [], "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                               "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": [],
                               "Anno_swiss": {"hit": [], "identity": []}, "Anno_nt": {"hit": [], "identity": []},
                               "Anno_nr": {"hit": [], "identity": []}, "Anno_neuropep": {"hit": [], "identity": []},
                               "Ortho": [], "ExSec": [], "Specificity": [],
                               "Whole_vs_Externa": [], "Whole_vs_Growing": [], "Whole_vs_Middle": [],
                               "Whole_vs_Terminal": [], "Growing_vs_Middle": [], "Growing_vs_Terminal": [],
                               "Middle_vs_Terminal": [], "Externa_vs_Growing": [], "Externa_vs_Middle": [],
                               "Externa_vs_Terminal": []}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".p")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"


def exsec_parsing(contig_dict, exsec_fasta, tag):
    exsec = SeqIO.parse(exsec_fasta, 'fasta')
    for seq in exsec:
        description = seq.id.strip().split(" ")
        exsec_id = description[0].split(".p")[0]
        contig_dict[exsec_id]["ExSec"].append(tag)


def orthogroups_parsing(contig_dict, orthogroups):
    head = orthogroups.readline()
    for line in orthogroups:
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        for protein_list in proteins:
            for protein in protein_list.split(","):
                if protein.split(".p")[0] in contig_dict.keys():
                    contig_dict[protein.split(".p")[0]]["Ortho"].append(group_ID)

    for contig, values in contig_dict.items():
        if len(values["Ortho"]) == 0:
            values["Ortho"].append("-")


def eggNOG_mapper_parsing(contig_dict, eggNOG):
    keys = ["EggNOG:Preferred_name", "EggNOG:GO_terms", "EggNOG:EC_number", "EggNOG:KEGG_KO", "EggNOG:KEGG_Pathway",
            "EggNOG:KEGG_Module", "EggNOG:KEGG_Reaction", "EggNOG:rclass", "EggNOG:BRITE", "EggNOG:KEGG_TC",
            "EggNOG:CAZy", "EggNOG:BiGG_Reaction", "EggNOG:OG", "EggNOG:COG_cat", "EggNOG:Description"]
    index_in_annotation = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20]
    for line in eggNOG:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            contig_ID, annotation = description[0].split(".p")[0], description[1:]
            for el in index_in_annotation:
                if len(annotation) >= el + 1 and len(annotation[el]) != 0:
                    contig_dict[contig_ID][keys[index_in_annotation.index(el)]].append(annotation[el])

    for contig_ID, values in contig_dict.items():
        for key in keys:
            if len(values[key]) == 0:
                values[key].append('-')


def BLAST_annotation(contig_dict, BLAST, key_tag):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        if key_tag == "swiss" or key_tag == "nr" or key_tag == "neuropep":
            ID, hit_name = description[0].strip().split(".p")[0], description[3]
        else:
            ID, hit_name = description[0].strip().split(" ")[0], description[3]

        if ID in contig_dict.keys():
            if hit_name == "no hits":
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append("No hit")
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append("No hit")
            else:
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["hit"].append(description[4])
                contig_dict[ID]["Anno_{key}".format(key=key_tag)]["identity"].append(description[9])

    if key_tag == "swiss" or key_tag == "nr" or key_tag == "neuropep":
        for contig in contig_dict.keys():
            if contig_dict[contig]["amino"] == "without ORF":
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["hit"].append("-")
                contig_dict[contig]["Anno_{key}".format(key=key_tag)]["identity"].append("-")


def domains_parsing(contig_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        ID, domain_name = description[0].strip().split(".p")[0], description[1]
        contig_dict[ID]["domains"].append(domain_name)

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["domains"]) == 0:
            contig_dict[contig]["domains"].append("-")


def jongeneel_specificity(contig_dict, file, tag):
    for line in file:
        contig = line.strip().split("\t")[0]
        contig_dict[contig]["Specificity"].append(tag)


def LFC_for_significant(contig_dict, table, tag):
    header = table.readline()
    for line in table:
        description = line.strip().split("\t")
        contig_ID, LFC = description[0], float(description[2])
        contig_dict[contig_ID][tag].append(LFC)

    for contig, values in contig_dict.items():
        if len(values[tag]) == 0:
            values[tag].append('-')


def write_output(contig_dict, out):
    with open("{out}.tsv".format(out=out), 'a') as output:
        output.write("### Contig_ID - ID of assembled supercontig\n"
                     "### Orthogroups - ID of the orthogroup constructed for analyzed Crustacea species\n"
                     "### NCBInt - best hit with NCBI nucleotide database\n"
                     "### NCBInr - best hit with NCBI non-redundant protein database\n"
                     "### SwissProt - best hit with SwissProt protein database\n"
                     "### NeuroPep - best hit with neuropeptide database\n"
                     "### NCBInt|NCBInr|SwissProt|NeuroPep_identity - best hit identity percentage\n"
                     "### Domain_arch - reconstructed protein architecture based on analysis "
                     " of comparison results with a PfamA database\n"
                     "### EggNOG:Preferred_name - Predicted gene name\n"
                     "### EggNOG:GO_terms - comma delimited list of predicted Gene Ontology terms\n"
                     "### EggNOG:EC_number - the Enzyme Commission number (EC number)\n"
                     "### EggNOG:KEGG_KO - comma delimited list of predicted KEGG KOs\n"
                     "### EggNOG:KEGG_Pathway - comma delimited list of predicted KEGG Pathways\n"
                     "### EggNOG:KEGG_Module - predicted KEGG Module\n"
                     "### EggNOG:KEGG_Reaction - predicted KEGG Reaction\n"
                     "### EggNOG:rclass - predicted KEGG Reaction Class\n"
                     "### EggNOG:BRITE - predicted KEGG BRITE\n"
                     "### EggNOG:KEGG_TC - predicted KEGG TC\n"
                     "### EggNOG:CAZy - best hit with Carbohydrate-Active enZYmes Database\n"
                     "### EggNOG:BiGG_Reaction - comma delimited list of predicted BiGG metabolic reactions\n"
                     "### EggNOG:OG - comma delimited list of matching eggNOG Orthologous Groups\n"
                     "### EggNOG:COG_cat - COG functional category inferred from best matching OG\n"
                     "### EggNOG:Description - eggNOG functional description inferred from best matching OG\n"
                     "### ExSec - classification of excretory/secretory proteins\n"
                     "### Specificity - expression specificity according to Jongeneel criterion\n"
                     "### LFC:Whole_vs_Externa|Growing|Middle|Terminal - expression Log2 Fold Change between "
                     "Whole body and Externa or Growing|Middle|Terminal part of stolon, respectively\n"
                     "### LFC:Growing_vs_Middle|Terminal and LFC:Middle_vs_Terminal - expression Log2 Fold Change "
                     "between different part of stolon\n"
                     "### LFC:Externa_vs_Growing|Middle|Terminal - expression Log2 Fold Change between Externa and "
                     "Growing|Middle|Terminal part of stolon, respectively\n"
                     "Contig_ID\tOrthogroup\tNCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                     "SwissProt\tSwissProt_identity\tNeuroPep\tNeuroPep_identity\tDomain_arch\tEggNOG:GO_terms\t"
                     "EggNOG:EC_number\t"
                     "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                     "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                     "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tExSec\tSpecificity\t"
                     "LFC:Whole_vs_Externa\t"
                     "LFC:Whole_vs_Growing\tLFC:Whole_vs_Middle\tLFC:Whole_vs_Terminal\tLFC:Growing_vs_Middle\t"
                     "LFC:Growing_vs_Terminal\tLFC:Middle_vs_Terminal\tLFC:Externa_vs_Growing\tLFC:Externa_vs_Middle\t"
                     "LFC:Externa_vs_Terminal\n")
        for contig, values in contig_dict.items():
            output.write("{id}\t{ortho}\t{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t{swiss}\t{swiss_identity}\t"
                         "{neuropep}\t{neuropep_identity}\t{domains}\t{go}\t{ec}\t{ko}\t{pathway}\t{module}\t"
                         "{reaction}\t{rclass}\t{brite}\t{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\t{exsec}\t"
                         "{specificity}\t{whole_vs_externa}\t{whole_vs_growing}\t{whole_vs_middle}\t"
                         "{whole_vs_terminal}\t{growing_vs_middle}\t{growing_vs_terminal}\t{middle_vs_terminal}\t"
                         "{externa_vs_growing}\t{externa_vs_middle}\t{externa_vs_terminal}\n".format(id=contig,
                          ortho=values["Ortho"][0], nt=values["Anno_nt"]["hit"][0],
                          nt_identity=values["Anno_nt"]["identity"][0],
                          nr=values["Anno_nr"]["hit"][0], nr_identity=values["Anno_nr"]["identity"][0],
                          swiss=values["Anno_swiss"]["hit"][0], swiss_identity=values["Anno_swiss"]["identity"][0],
                          neuropep=values["Anno_neuropep"]["hit"][0],
                          neuropep_identity=values["Anno_neuropep"]["identity"][0],
                          domains="|".join(values["domains"]), go=values["EggNOG:GO_terms"][0],
                          ec=values["EggNOG:EC_number"][0], ko=values["EggNOG:KEGG_KO"][0],
                          pathway=values["EggNOG:KEGG_Pathway"][0], module=values["EggNOG:KEGG_Module"][0],
                          reaction=values["EggNOG:KEGG_Reaction"][0], rclass=values["EggNOG:rclass"][0],
                          brite=values["EggNOG:BRITE"][0], tc=values["EggNOG:KEGG_TC"][0],
                          cazy=values["EggNOG:CAZy"][0],
                          bigg=values["EggNOG:BiGG_Reaction"][0], og=values["EggNOG:OG"][0],
                          cog=values["EggNOG:COG_cat"][0], description=values["EggNOG:Description"][0],
                          exsec=values["ExSec"][0], specificity=values["Specificity"][0],
                          whole_vs_externa=values["Whole_vs_Externa"][0],
                          whole_vs_growing=values["Whole_vs_Growing"][0],
                          whole_vs_middle=values["Whole_vs_Middle"][0],
                          whole_vs_terminal=values["Whole_vs_Terminal"][0],
                          growing_vs_middle=values["Growing_vs_Middle"][0],
                          growing_vs_terminal=values["Growing_vs_Terminal"][0],
                          middle_vs_terminal=values["Middle_vs_Terminal"][0],
                          externa_vs_growing=values["Externa_vs_Growing"][0],
                          externa_vs_middle=values["Externa_vs_Middle"][0],
                          externa_vs_terminal=values["Externa_vs_Terminal"][0]))


if __name__ == "__main__":
        contig_dict = {}
        print("***** Fasta files parsing *****")
        nucl_parsing(contig_dict, args.nucl)
        amino_parsing(contig_dict, args.amino)
        exsec_parsing(contig_dict, args.class_exsec, "Classical")
        exsec_parsing(contig_dict, args.nonclass_exsec, "Non-classical")
        for contig, values in contig_dict.items():
            if len(values["ExSec"]) == 0:
                values["ExSec"].append('-')
        print("***** Orthogpoups parsing *****")
        orthogroups_parsing(contig_dict, args.ortho)
        print("***** EggNOG-mapper output parsing *****")
        eggNOG_mapper_parsing(contig_dict, args.eggnog)
        print("***** BLAST results parsing *****")
        BLAST_annotation(contig_dict, args.blast_swiss, "swiss")
        BLAST_annotation(contig_dict, args.blast_nt, "nt")
        BLAST_annotation(contig_dict, args.blast_nr, "nr")
        BLAST_annotation(contig_dict, args.neuropep, "neuropep")
        print("***** Domain architecture parsing *****")
        domains_parsing(contig_dict, args.domains)
        print("***** Specificity searching output files parsing *****")
        jongeneel_specificity(contig_dict, args.externa_jongeneel, "Externa-specific")
        jongeneel_specificity(contig_dict, args.growing_jongeneel, "Growing_stolon-specific")
        jongeneel_specificity(contig_dict, args.middle_jongeneel, "Middle_stolon_part-specific")
        jongeneel_specificity(contig_dict, args.terminal_jongeneel, "Terminal_stolon_part-specific")
        jongeneel_specificity(contig_dict, args.whole_jongeneel, "Whole_body-specific")
        for contig, values in contig_dict.items():
            if len(values["Specificity"]) == 0:
                values["Specificity"].append('-')
        print("***** DESeq2-contrast (significant) output parsing *****")
        LFC_for_significant(contig_dict, args.whole_vs_externa_sign, "Whole_vs_Externa")
        LFC_for_significant(contig_dict, args.whole_vs_growing_sign, "Whole_vs_Growing")
        LFC_for_significant(contig_dict, args.whole_vs_middle_sign, "Whole_vs_Middle")
        LFC_for_significant(contig_dict, args.whole_vs_terminal_sign, "Whole_vs_Terminal")
        LFC_for_significant(contig_dict, args.growing_vs_middle_sign, "Growing_vs_Middle")
        LFC_for_significant(contig_dict, args.growing_vs_terminal_sign, "Growing_vs_Terminal")
        LFC_for_significant(contig_dict, args.middle_vs_terminal_sign, "Middle_vs_Terminal")
        LFC_for_significant(contig_dict, args.externa_vs_growing_sign, "Externa_vs_Growing")
        LFC_for_significant(contig_dict, args.externa_vs_middle_sign, "Externa_vs_Middle")
        LFC_for_significant(contig_dict, args.externa_vs_terminal_sign, "Externa_vs_Terminal")
        print("***** Output file creating *****")
        write_output(contig_dict, args.out)
