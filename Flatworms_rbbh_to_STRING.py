try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()


parser = argparse.ArgumentParser()
parser.add_argument('--csin', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Csinensis and STRINGdb")
parser.add_argument('--psim', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Psimillimum and STRINGdb")
parser.add_argument('--fhep', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Fhepatica and STRINGdb")
parser.add_argument('--fgig', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Fgigantica and STRINGdb")
parser.add_argument('--ofel', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Ofelineus and STRINGdb")
parser.add_argument('--oviv', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Oviverrini and STRINGdb")
parser.add_argument('--shae', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Shaematobium and STRINGdb")
parser.add_argument('--sjap', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Sjaponicum and STRINGdb")
parser.add_argument('--sman', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Smansoni and STRINGdb")
parser.add_argument('--treg', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Tregenti and STRINGdb")
parser.add_argument('--tszi', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Tszidati and STRINGdb")
parser.add_argument('--smed', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Smediterranea and STRINGdb")
parser.add_argument('--mlig', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Mlignano and STRINGdb")
parser.add_argument('--pvit', type=argparse.FileType('r'), required=True,
                    help="Table with reciprocal best BLAST hits between Pvitattus and STRINGdb")
parser.add_argument('--string_tag', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


def rbbh_parsing(rbbh_tab, rbbh_dict, species_tag):
    header = rbbh_tab.readline()
    for line in rbbh_tab:
        description = line.strip().split("\t")
        rbbh_ID, species_seq, string_seq = description[0], description[1], description[2]
        if string_seq not in rbbh_dict.keys():
            rbbh_dict[string_seq] = {"csin": [], "fgig": [], "fhep": [], "ofel": [], "oviv": [], "psim": [],
                                     "sjap": [], "sman": [], "shae": [], "treg": [], "tszi": [],
                                     "smed": [], "mlig": [], "pvit": []}
            rbbh_dict[string_seq][species_tag].append(species_seq)
        else:
            rbbh_dict[string_seq][species_tag].append(species_seq)


def output_writing(rbbh_dict, string_tag, output):
    with open("{output}.all_rbbh.tsv".format(output=output), 'a') as all_rbbh:
        all_rbbh.write("{string}\tCsinensis\tFgigantica\tFhepatica\tOfelineus\tOviverrini\tPsimillimum\t"
                       "Shaematobium\tSjaponicum\tSmansoni\tTregenti\tTszidati\tSmediterranea\tMlignano\t"
                       "Pvittatus\n".format(string=string_tag))
        for string_seq, values in rbbh_dict.items():
            all_rbbh.write("{string}\t{csin}\t{fgig}\t{fhep}\t{ofel}\t{oviv}\t{psim}\t"
                           "{shae}\t{sjap}\t{sman}\t{treg}\t{tszi}\t{smed}\t{mlig}\t{pvit}\n".format(
                            string=string_seq, csin=values["csin"][0], fgig=values["fgig"][0], fhep=values["fhep"][0],
                            ofel=values["ofel"][0], oviv=values["oviv"][0], psim=values["psim"][0],
                            shae=values["shae"][0], sjap=values["sjap"][0], sman=values["sman"][0],
                            treg=values["treg"][0], tszi=values["tszi"][0], smed=values["smed"][0],
                            mlig=values["mlig"][0], pvit=values["pvit"][0]))

    with open("{output}.rbbh_common_for_all_flatworms.tsv".format(output=output), 'a') as flatworm_rbbh:
        flatworm_rbbh.write("{string}\tCsinensis\tFgigantica\tFhepatica\tOfelineus\tOviverrini\tPsimillimum\t"
                            "Shaematobium\tSjaponicum\tSmansoni\tTregenti\tTszidati\tSmediterranea\tMlignano\t"
                            "Pvittatus\n".format(string=string_tag))
        for string_seq, values in rbbh_dict.items():
            hits = [values[species][0] for species in ["csin", "fgig", "fhep", "ofel", "oviv", "psim", "sjap",
                                                       "sman", "shae", "treg", "tszi", "smed", "mlig", "pvit"]]
            if hits.count("-") == 0:
                flatworm_rbbh.write("{string}\t{csin}\t{fgig}\t{fhep}\t{ofel}\t{oviv}\t{psim}\t"
                                    "{shae}\t{sjap}\t{sman}\t{treg}\t{tszi}\t{smed}\t{mlig}\t{pvit}\n".format(
                                     string=string_seq, csin=values["csin"][0], fgig=values["fgig"][0],
                                     fhep=values["fhep"][0], ofel=values["ofel"][0], oviv=values["oviv"][0],
                                     psim=values["psim"][0], shae=values["shae"][0], sjap=values["sjap"][0],
                                     sman=values["sman"][0], treg=values["treg"][0], tszi=values["tszi"][0],
                                     smed=values["smed"][0], mlig=values["mlig"][0], pvit=values["pvit"][0]))

    with open("{output}.rbbh_common_for_trematodes.tsv".format(output=output), 'a') as trematoda_rbbh:
        trematoda_rbbh.write("{string}\tCsinensis\tFgigantica\tFhepatica\tOfelineus\tOviverrini\tPsimillimum\t"
                             "Shaematobium\tSjaponicum\tSmansoni\tTregenti\tTszidati\n".format(string=string_tag))
        for string_seq, values in rbbh_dict.items():
            hits = [values[species][0] for species in ["csin", "fgig", "fhep", "ofel", "oviv", "psim", "sjap",
                                                       "sman", "shae", "treg", "tszi"]]
            if hits.count("-") == 0:
                trematoda_rbbh.write("{string}\t{csin}\t{fgig}\t{fhep}\t{ofel}\t{oviv}\t{psim}\t"
                                     "{shae}\t{sjap}\t{sman}\t{treg}\t{tszi}\n".format(
                                      string=string_seq, csin=values["csin"][0], fgig=values["fgig"][0],
                                      fhep=values["fhep"][0], ofel=values["ofel"][0], oviv=values["oviv"][0],
                                      psim=values["psim"][0], shae=values["shae"][0], sjap=values["sjap"][0],
                                      sman=values["sman"][0], treg=values["treg"][0], tszi=values["tszi"][0]))

    with open("{output}.rbbh_common_for_free-living.tsv".format(output=output), 'a') as freeliving_rbbh:
        freeliving_rbbh.write("{string}\tSmediterranea\tMlignano\tPvittatus\n".format(string=string_tag))
        for string_seq, values in rbbh_dict.items():
            hits = [values[species][0] for species in ["smed", "mlig", "pvit"]]
            if hits.count("-") == 0:
                freeliving_rbbh.write("{string}\t{smed}\t{mlig}\t{pvit}\n".format(string=string_seq,
                                      smed=values["smed"][0], mlig=values["mlig"][0], pvit=values["pvit"][0]))


if __name__ == "__main__":
    rbbh_dict = {}
    print("***** Input files parsing *****")
    rbbh_parsing(args.csin, rbbh_dict, "csin")
    rbbh_parsing(args.psim, rbbh_dict, "psim")
    rbbh_parsing(args.fgig, rbbh_dict, "fgig")
    rbbh_parsing(args.fhep, rbbh_dict, "fhep")
    rbbh_parsing(args.ofel, rbbh_dict, "ofel")
    rbbh_parsing(args.oviv, rbbh_dict, "oviv")
    rbbh_parsing(args.shae, rbbh_dict, "shae")
    rbbh_parsing(args.sjap, rbbh_dict, "sjap")
    rbbh_parsing(args.sman, rbbh_dict, "sman")
    rbbh_parsing(args.treg, rbbh_dict, "treg")
    rbbh_parsing(args.tszi, rbbh_dict, "tszi")
    rbbh_parsing(args.smed, rbbh_dict, "smed")
    rbbh_parsing(args.mlig, rbbh_dict, "mlig")
    rbbh_parsing(args.pvit, rbbh_dict, "pvit")

    for string_seq, values in rbbh_dict.items():
        for species_tag, value in values.items():
            if len(value) == 0:
                value.append("-")

    print("***** Output files writing *****")
    output_writing(rbbh_dict, args.string_tag, args.output)


