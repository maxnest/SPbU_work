try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=str)
args = parser.parse_args()


def domains_parsing(domains_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        protein_ID, domain_name, domain_ID, start, end = description[0], description[1], description[2], \
                                                         description[3], description[4]
        domain_info = "{id}:{name};{start}:{end}".format(name=domain_name, id=domain_ID, start=start, end=end)
        if protein_ID not in domains_dict.keys():
            domains_dict[protein_ID] = []
        domains_dict[protein_ID].append(domain_info)


def output_writing(domains_dict, out):
    with open("{out}.tsv".format(out=out), 'a') as output_file:
        output_file.write("Protein_IDs\tDomain_arch\n")
        for protein_ID, domains in domains_dict.items():
            output_file.write("{id}\t{domain_arch}\n".format(id=protein_ID, domain_arch="|".join(domains)))


if __name__ == "__main__":
    domains_dict = {}
    print("***** Table parsing *****")
    domains_parsing(domains_dict, args.domains)
    print("***** Output file writing *****")
    output_writing(domains_dict, args.out)