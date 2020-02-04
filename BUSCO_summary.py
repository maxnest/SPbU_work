try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--busco_all_trinity', type=argparse.FileType('r'), required=True)
parser.add_argument('--busco_good_trinity', type=argparse.FileType('r'), required=True)
parser.add_argument('--busco_all_cdhit_sense', type=argparse.FileType('r'), required=True)
parser.add_argument('--busco_good_cdhit_sense', type=argparse.FileType('r'), required=True)
parser.add_argument('--busco_good_cdhit_sense_pep', type=argparse.FileType('r'), required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def busco_full_tab_parsing(busco_dict, busco_tab, tag):
    for line in busco_tab:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            ID, status = description[0], description[1]
            if ID not in busco_dict.keys():
                keys = ["status", "sequences", "score"]
                busco_dict[ID] = {"all_Trinity": {key: [] for key in keys},
                                  # "all_Trinity" - all assembled sequences (Trinity, before clusterization)
                                  "good_Trinity": {key: [] for key in keys},
                                  # "good_Trinity" - only sequences, classified as "good" by TransRate
                                  "all_cdhit_sense": {key: [] for key in keys},
                                  # "all_cdhit_sense" - all sequences selected by CDHIT (only +/+ mode)
                                  "good_cdhit_sense": {key: [] for key in keys},
                                  # "good_cdhit_sense" - only sequences, classified as "good" by TransRate
                                  "good_cdhit_sense_pep": {key: [] for key in keys}
                                  # "good_cdhit_sense_pep" - proteins, predicted on the clustered and selected contigs
                }

            if status != "Missing":
                busco_dict[ID][tag]["status"].append(status)
                busco_dict[ID][tag]["sequences"].append(description[2])
                busco_dict[ID][tag]["score"].append(float(description[3]))
            else:
                busco_dict[ID][tag]["status"].append(status)
                busco_dict[ID][tag]["score"].append(0)


def output_writing(busco_dict, output):
    with open("{output}.tsv".format(output=output), 'a') as output_file:
        output_file.write("BUSCO ID\tTrinity_all_seq:status\tTrinity_only_good:status\t"
                          "CDHIT_all_seq:status\tCDHIT_only_good:status\tCDHIT_pep:status\t"
                          "Trinity_all:min_and_max_scores\tTrinity_only_good:min_and_max_scores\t"
                          "CDHIT_all_seq:min_and_max_scores\tCDHIT_only_good:min_and_max_scores\t"
                          "CDHIT_pep:min_and_max_scores\n")
        for id, values in busco_dict.items():
            output_file.write("{id}\t{all_status}\t{good_status}\t{cdhit_all_status}\t"
                              "{cdhit_good_status}\t{cdhit_pep_status}\t{all_scores}\t{good_scores}\t"
                              "{cdhit_all_scores}\t{cdhit_good_scores}\t{cdhit_pep_scores}\n".format(
                               id=id, all_status=list(set(values["all_Trinity"]["status"]))[0],
                               good_status=list(set(values["good_Trinity"]["status"]))[0],
                               cdhit_all_status=list(set(values["all_cdhit_sense"]["status"]))[0],
                               cdhit_good_status=list(set(values["good_cdhit_sense"]["status"]))[0],
                               cdhit_pep_status=list(set(values["good_cdhit_sense_pep"]["status"]))[0],
                               all_scores="{min}|{max}".format(min=np.min(values["all_Trinity"]["score"]),
                                                               max=np.max(values["all_Trinity"]["score"])),
                               good_scores="{min}|{max}".format(min=np.min(values["good_Trinity"]["score"]),
                                                                max=np.max(values["good_Trinity"]["score"])),
                               cdhit_all_scores="{min}|{max}".format(min=np.min(values["all_cdhit_sense"]["score"]),
                                                                     max=np.max(values["all_cdhit_sense"]["score"])),
                               cdhit_good_scores="{min}|{max}".format(min=np.min(values["good_cdhit_sense"]["score"]),
                                                                      max=np.max(values["good_cdhit_sense"]["score"])),
                               cdhit_pep_scores="{min}|{max}".format(
                                min=np.min(values["good_cdhit_sense_pep"]["score"]),
                                max=np.max(values["good_cdhit_sense_pep"]["score"]))))


if __name__ == "__main__":
    busco_dict = {}
    print("***** Input files parsing *****")
    busco_full_tab_parsing(busco_dict, args.busco_all_trinity, "all_Trinity")
    busco_full_tab_parsing(busco_dict, args.busco_good_trinity, "good_Trinity")
    busco_full_tab_parsing(busco_dict, args.busco_all_cdhit_sense, "all_cdhit_sense")
    busco_full_tab_parsing(busco_dict, args.busco_good_cdhit_sense, "good_cdhit_sense")
    busco_full_tab_parsing(busco_dict, args.busco_good_cdhit_sense_pep, "good_cdhit_sense_pep")
    print("***** Output file creating *****")
    output_writing(busco_dict, args.output)
