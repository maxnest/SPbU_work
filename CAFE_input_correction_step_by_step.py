try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

import os
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=argparse.FileType('r'), required=True)
parser.add_argument('--sh_file', type=argparse.FileType('r'), required=True)
args = parser.parse_args()


def read_input(input, input_dict):
    head = input.readline()
    for line in input:
        description = line.strip().split("\t")
        cluster, group, csin, fhep, mlig, oviv, shae, sjap, sman, smed = description[0], description[1], \
                                                                         int(description[2]), int(description[3]), \
                                                                         int(description[4]), int(description[5]), \
                                                                         int(description[6]), int(description[7]), \
                                                                         int(description[8]), int(description[9])

        input_dict[group] = {"Cluster": cluster, "Csin": csin, "Fhep": fhep, "Mlig": mlig, "Oviv": oviv, "Shae": shae,
                             "Sjap": sjap, "Sman": sman, "Smed": smed}


def read_sh(sh, sh_list):
    for line in sh:
        sh_list.append(line)


def sh_alter(sh_list, count):
    with open("test_{count}.sh".format(count=count), 'a') as new_sh:
        new_sh.write(sh_list[0])
        new_sh.write("load -i Test_matrix_{count}.tsv -t 15 -l test_{count}.log\n".format(count=count))
        new_sh.write("{second_line}\n{third_line}\n".format(second_line=sh_list[2], third_line=sh_list[3]))
        new_sh.write("report cafe_reports/test_{count}.report\n".format(count=count))
    new_sh.close()


def writing_new_matrix(count, warnings, input_dict):
    if len(warnings) != len(input_dict.keys()):
        with open("Test_matrix_{count}.tsv".format(count=count), 'a') as new_matrix:
            new_matrix.write("Description\tID\tCsin\tFhep\tMlig\tOviv\tShae\tSjap\tSman\tSmed\n")
            for group, values in input_dict.items():
                if group not in warnings:
                    new_matrix.write("{cluster}\t{group}\t{csin}\t{fhep}\t{mlig}\t{oviv}\t"
                                     "{shae}\t{sjap}\t{sman}\t{smed}\n".format(cluster=values["Cluster"], group=group,
                                                                               csin=values["Csin"], fhep=values["Fhep"],
                                                                               mlig=values["Mlig"], oviv=values["Oviv"],
                                                                               shae=values["Shae"], sjap=values["Sjap"],
                                                                               sman=values["Sman"], smed=values["Smed"]))


if __name__ == "__main__":
    cwd = os.getcwd()
    count = 0
    input_dict, warnings, sh_list = {}, [], []
    print("***** Files reading *****")
    read_input(args.input, input_dict)
    read_sh(args.sh_file, sh_list)
    print("***** Cafe running *****")
    while len(warnings) >= count:
        writing_new_matrix(count, warnings, input_dict)
        sh_alter(sh_list, count)
        cafe_run = "cafe {sh}".format(sh="{cwd}/test_{count}.sh".format(cwd=cwd, count=count))
        process = subprocess.Popen(cafe_run.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            for line in process.stderr:
                if line.decode().split()[-3] not in warnings:
                    warnings.append(line.decode().split()[-3])
            print("***** {warning} was excluded *****".format(warning=warnings[-1]))
        except:
            process.kill()
            output, error = process.communicate()
        count += 1

    with open("warnings.txt", 'a') as warnings_output:
        for el in warnings:
            warnings_output.write("{el}\n".format(el=el))