try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_swiss', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nt', type=argparse.FileType('r'), required=True)
parser.add_argument('--blast_nr', type=argparse.FileType('r'), required=True)
parser.add_argument('--eggnog', type=argparse.FileType('r'), required=True)
parser.add_argument('--sites_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_sign', type=argparse.FileType('r'), required=True)
parser.add_argument('--head_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--tail_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--sites_head_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--sites_tail_clusters', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id.split(" ")[0]] = {"ID": [seq.id.split(" ")[0]], "nucl": seq.seq, "amino": "",
                                             "domains": [], "Anno_swiss": {"hit": [], "identity": []},
                                             "Anno_nt": {"hit": [], "identity": []},
                                             "Anno_nr": {"hit": [], "identity": []},
                                             "head_significant": [], "tail_significant": [], "sites_significant": [],
                                             "head_cluster": [], "tail_cluster": [], "sites_head_cluster": [],
                                             "sites_tail_cluster": [],
                                             "EggNOG:Preferred_name": [], "EggNOG:GO_terms": [], "EggNOG:EC_number": [],
                                             "EggNOG:KEGG_KO": [], "EggNOG:KEGG_Pathway": [], "EggNOG:KEGG_Module": [],
                                             "EggNOG:KEGG_Reaction": [], "EggNOG:rclass": [], "EggNOG:BRITE": [],
                                             "EggNOG:KEGG_TC": [], "EggNOG:CAZy": [], "EggNOG:BiGG_Reaction": [],
                                             "EggNOG:OG": [], "EggNOG:COG_cat": [], "EggNOG:Description": []}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".p")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"


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
        if key_tag == "swiss" or key_tag == "nr":
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

    if key_tag == "swiss" or key_tag == "nr":
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


def clusters(contig_dict, table, dict, tag):
    header = table.readline().strip().split("\t")

    for el in header:
        dict[el] = []

    for line in table:
        genes = line.strip().split("\t")
        for gene in genes:
            if len(gene) != 0:
                dict[header[genes.index(gene)]].append(gene)

    for contig in contig_dict.keys():
        for cluster, genes in dict.items():
            if contig in genes:
                contig_dict[contig]["{tag}_cluster".format(tag=tag)].append(cluster)

    for contig, values in contig_dict.items():
        if len(values["{tag}_cluster".format(tag=tag)]) == 0:
            values["{tag}_cluster".format(tag=tag)].append("-")


def significant(contig_dict, table, tag):
    for line in table:
        contig = line.strip().split(",")[0]
        contig_dict[contig]["{tag}_significant".format(tag=tag)].append("*")

    for contig, values in contig_dict.items():
        if len(values["{tag}_significant".format(tag=tag)]) == 0:
            values["{tag}_significant".format(tag=tag)].append("-")


def write_output_files(contig_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("### Contig_ID - ID of assembled supercontig\n"
                          "### NCBInt - best hit with NCBI nucleotide database\n"
                          "### NCBInr - best hit with NCBI non-redundant protein database\n"
                          "### SwissProt - best hit with SwissProt protein database\n"
                          "### NCBInt|NCBInr|SwissProt_identity - best hit identity percentage\n"
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
                          "### Sites_significant - the sequence expression differs significantly between "
                          "the two sites (Head vs Tail)\n"
                          "### Sites_Head_cluster - Cluster ID of site-significant co-expressed genes "
                          "reconstructed based on the values of gene activity in the Head\n"
                          "### Sites_Tail_cluster - Cluster ID of site-significant co-expressed genes "
                          "reconstructed based on the values of gene activity in the Tail\n"
                          "### Head|Tail_significant - the sequence expression differs significantly between "
                          "the time points in Head or Tail, respectively\n"
                          "### Head_cluster - Cluster ID of Head-significant co-expressed genes \n"
                          "### Tail_cluster - Cluster ID of Tail-significant co-expressed genes\n")
        output_file.write("Contig_ID\tEggNOG:Preferred_name\tNCBInt\tNCBInt_identity\tNCBInr\tNCBInr_identity\t"
                          "SwissProt\tSwissProt_identity\tDomains_arch\tEggNOG:GO_terms\tEggNOG:EC_number\t"
                          "EggNOG:KEGG_KO\tEggNOG:KEGG_Pathway\tEggNOG:KEGG_Module\tEggNOG:KEGG_Reaction\t"
                          "EggNOG:rclass\tEggNOG:BRITE\tEggNOG:KEGG_TC\tEggNOG:CAZy\tEggNOG:BiGG_Reaction\t"
                          "EggNOG:OG\tEggNOG:COG_cat\tEggNOG:Description\tSites_significant\tSites_Head_cluster\t"
                          "Sites_Tail_cluster\tHead_significant\tHead_cluster\tTail_significant\tTail_cluster\n")
        for contig, values in contig_dict.items():
            output_file.write("{id}\t{name}\t{nt}\t{nt_identity}\t{nr}\t{nr_identity}\t{swiss}\t{swiss_identity}\t"
                              "{domains}\t{go}\t{ec}\t{ko}\t{pathway}\t{module}\t{reaction}\t{rclass}\t{brite}\t"
                              "{tc}\t{cazy}\t{bigg}\t{og}\t{cog}\t{description}\t{sites_sign}\t"
                              "{sites_head_cluster}\t{sites_tail_cluster}\t{head_sign}\t{head_cluster}\t{tail_sign}\t"
                              "{tail_cluster}\n".format(
                               id=contig, name=values["EggNOG:Preferred_name"][0], nt=values["Anno_nt"]["hit"][0],
                               nt_identity=values["Anno_nt"]["identity"][0], nr=values["Anno_nr"]["hit"][0],
                               nr_identity=values["Anno_nr"]["identity"][0], swiss=values["Anno_swiss"]["hit"][0],
                               swiss_identity=values["Anno_swiss"]["identity"][0], domains="|".join(values["domains"]),
                               go=values["EggNOG:GO_terms"][0], ec=values["EggNOG:EC_number"][0],
                               ko=values["EggNOG:KEGG_KO"][0], pathway=values["EggNOG:KEGG_Pathway"][0],
                               module=values["EggNOG:KEGG_Module"][0], reaction=values["EggNOG:KEGG_Reaction"][0],
                               rclass=values["EggNOG:rclass"][0], brite=values["EggNOG:BRITE"][0],
                               tc=values["EggNOG:KEGG_TC"][0], cazy=values["EggNOG:CAZy"][0],
                               bigg=values["EggNOG:BiGG_Reaction"][0], og=values["EggNOG:OG"][0],
                               cog=values["EggNOG:COG_cat"][0], description=values["EggNOG:Description"][0],
                               sites_sign=values["sites_significant"][0],
                               sites_head_cluster=values["sites_head_cluster"][0],
                               sites_tail_cluster=values["sites_tail_cluster"][0],
                               head_sign=values["head_significant"][0], head_cluster=values["head_cluster"][0],
                               tail_sign=values["tail_significant"][0], tail_cluster=values["tail_cluster"][0]))


if __name__ == "__main__":
    contig_dict = {}
    head_clusters, tail_clusters, sites_head_clusters, sites_tail_clusters = {}, {}, {}, {}
    print("***** Fasta files parsing *****")
    nucl_parsing(contig_dict, args.nucl)
    amino_parsing(contig_dict, args.amino)
    print("***** eggNOG-mapper output parsing *****")
    eggNOG_mapper_parsing(contig_dict, args.eggnog)
    print("***** BLAST results parsing *****")
    BLAST_annotation(contig_dict, args.blast_swiss, "swiss")
    BLAST_annotation(contig_dict, args.blast_nt, "nt")
    BLAST_annotation(contig_dict, args.blast_nr, "nr")
    print("***** Domain architecture parsing *****")
    domains_parsing(contig_dict, args.domains)
    print("***** Cluster parsing *****")
    clusters(contig_dict, args.head_clusters, head_clusters, "head")
    clusters(contig_dict, args.tail_clusters, tail_clusters, "tail")
    clusters(contig_dict, args.sites_head_clusters, sites_head_clusters, "sites_head")
    clusters(contig_dict, args.sites_tail_clusters, sites_tail_clusters, "sites_tail")
    print("***** Tables with IDs of significant sequences parsing *****")
    significant(contig_dict, args.sites_sign, "sites")
    significant(contig_dict, args.head_sign, "head")
    significant(contig_dict, args.tail_sign, "tail")
    print("***** Output file creating *****")
    write_output_files(contig_dict, args.output)
