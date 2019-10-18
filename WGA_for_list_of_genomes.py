try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

from subprocess import call

parser = argparse.ArgumentParser()
parser.add_argument('--genomes', type=argparse.FileType('r'), required=True,
                    help="Text file with names of fasta-files with genomes one per line")
parser.add_argument('--threads', type=int, required=True)
args = parser.parse_args()


def wga_for_pair(reference, query, threads):
    ref_tag, query_tag = reference.split("/")[-1].split(".")[0], query.split("/")[-1].split(".")[0]
    call(["mkdir", "{ref_tag}_vs_{query_tag}".format(ref_tag=ref_tag, query_tag=query_tag)])
    call("cd {ref_tag}_vs_{query_tag}".format(ref_tag=ref_tag, query_tag=query_tag), shell=True)
    call(["nohup", "nucmer", "-p", "{prefix}".format(prefix="{ref_tag}_vs_{query_tag}".format(
        ref_tag=ref_tag, query_tag=query_tag)), "-t", "{threads}".format(threads=threads),
          "{ref}".format(ref=reference), "{query}".format(query=query)])
    call(["nohup", "dnadiff", "-p", "{prefix}".format(prefix="{ref_tag}_vs_{query_tag}".format(
        ref_tag=ref_tag, query_tag=query_tag)), "-d", "{prefix}.delta".format(prefix="{ref_tag}_vs_{query_tag}".format(
        ref_tag=ref_tag, query_tag=query_tag))])
    call("cd ..", shell=True)


if __name__ == "__main__":
    genomes_list = []
    for line in args.genomes:
        genomes_list.append(line.strip())
    for ref_genome in genomes_list:
        ref_tag = ref_genome.split("/")[-1].split(".")[0]
        call(["mkdir", "{ref_tag}".format(ref_tag=ref_tag)])
        call("cd {ref_tag}/".format(ref_tag=ref_tag), shell=True)
        for query_genome in genomes_list:
            if ref_genome != query_genome:
                print("WGA for {ref} and {query}".format(ref=ref_genome, query=query_genome))
                wga_for_pair(ref_genome, query_genome, args.threads)
        call("cd ..", shell=True)
