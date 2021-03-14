try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import collections

parser = argparse.ArgumentParser()
parser.add_argument('--leapfrog_merged', type=argparse.FileType('r'), required=True,
                    help="The tsv file with merged leapfrog results (after using different species as the bridge)")
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def leapfrog_merged_parsing(leapfrog_merged, leapfrog_merged_dict):
    header = leapfrog_merged.readline()
    for line in leapfrog_merged:
        description = line.strip().split("\t")
        query, mid, subject, query_vs_mid, mid_vs_subject = description[0], description[1], description[2], \
                                                            float(description[3]), float(description[4])
        if query not in leapfrog_merged_dict:
            leapfrog_merged_dict[query] = {}

        leapfrog_merged_dict[query]["{mid}_vs_{subject}".format(mid=mid, subject=subject)] = {
            "mid_vs_subject": mid_vs_subject, "query_vs_mid": query_vs_mid}


def bbh_searching(leapfrog_merged_dict, leapfrog_bbh_dict):
    for query, hits in leapfrog_merged_dict.items():
        values = [value["mid_vs_subject"] for value in hits.values()]
        for hit, value in hits.items():
            if value["mid_vs_subject"] == min(values):
                leapfrog_bbh_dict[query] = {"mid": hit.split("_vs_")[0], "subject": hit.split("_vs_")[1],
                                            "mid_vs_subject": value["mid_vs_subject"],
                                            "query_vs_mid": value["query_vs_mid"]}
                break


def output_writing(output, leapfrog_bbh_dict):
    with open("{output}.leapfrog_merged_bbh.tsv".format(output=output), 'a') as output_file:
        output_file.write("Query\tMid\tSubject\tQuery.v.Mid.E.Val\tMid.v.Subject.E.Val\n")
        for query, value in leapfrog_bbh_dict.items():
            output_file.write("{query}\t{mid}\t{subject}\t{query_vs_mid}\t{mid_vs_subject}\n".format(
                query=query, mid=value["mid"], subject=value["subject"], query_vs_mid=value["query_vs_mid"],
                mid_vs_subject=value["mid_vs_subject"]))


if __name__ == "__main__":
    leapfrog_merged_dict, leapfrog_bbh_dict = {}, {}
    leapfrog_merged_parsing(args.leapfrog_merged, leapfrog_merged_dict)
    bbh_searching(leapfrog_merged_dict, leapfrog_bbh_dict)
    output_writing(args.output, leapfrog_bbh_dict)


