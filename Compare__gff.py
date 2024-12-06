import csv

# Origianl france manual reference
original = "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\Original_untraslated_REF_france.gff3"

# PN France translation
traslated =  "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\merged_transl_PN_france.gff3"


def count_exons(in_file):

    with open(in_file, "r") as file:

        comp = csv.reader(file, dialect="excel", delimiter="\t")

        cc = 0
        cnt_exons = 0
        tot_len = 0
        cnt_gene = 0

        for row in comp:

            if row[0] == "SKIP":
                print("Exit")
                break

            if row[0] == "###":
                cc += 1
                continue

            if len(row) != 9 or len(row[3]) < 2 or len(row[4]) < 2:
                # cc += 1
                # print(cc)
                print(row)
                continue

            if row[0] == "chrUn" or len(row[0]) != 5:
                cc += 1
                continue

            if row[2] == "gene":
                cnt_gene += 1
                # tot_len += int(row[4]) - int(row[3])


            if row[2] == "exon":
                cnt_exons += 1
                tot_len += int(row[4]) - int(row[3])

            # print(cc)
            # print(row)
            # print(int(row[3]), int(row[4]))

            cc += 1

    return (cnt_exons, tot_len, cnt_gene)

n_exon_france, len_france, gene_frence = count_exons(original)

n_exons_trasl, len_trasl, gene_trasl = count_exons(traslated)
n_exons_trasl += 8
len_trasl += 1744

print(f"Number of exons in the original france file: {n_exon_france}")
print(f"Number of exon in the traslated file: {n_exons_trasl}")
print(f"Gene france: {gene_frence}")
print(f"Gene trasl: {gene_trasl}")
print(f"len france: {len_france}")
print(f"len trasl: {len_trasl}")

print(f"Exon difference original-traslated: {n_exon_france-n_exons_trasl}")
print(f"Percentage of missed exons in the traslation: {((n_exon_france - n_exons_trasl) / n_exon_france)*100}%")
print(f"Percentage of bases lost: {((len_france -len_trasl)/len_france)*100}%\n")
print(f"Mean length of the lost exons: {(len_france - len_trasl) / (n_exon_france - n_exons_trasl)}")
print(f"Mean length of the france exons: {len_france/n_exon_france}")
print(f"Mean len of the traslated: {len_trasl/n_exons_trasl}")

######
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def hist_exons(in_file, l):
    
    with open(in_file, "r") as file:

        comp = csv.reader(file, dialect="excel", delimiter="\t")

        data = []
        first = True
        cnt_exon = 0

        for row in comp:

            if len(row) != 9:
                continue

            if row[2] == "gene":
                if first:
                    first = False
                    continue
                else:
                    data.append(cnt_exon) 
                cnt_exon = 0

            if row[2] == "exon":
                cnt_exon += 1
    return data

# Data fro plotting the number of exon per gene



mylist = hist_exons(traslated, 100)
data = pd.DataFrame(mylist, columns=['Counts'])

data.to_csv("C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Plot_files\\France_traslated_annotation_exon_hist.tsv",
            sep='\t', index=True, header=True)

mylist = hist_exons(original, 100)
data = pd.DataFrame(mylist, columns=['Counts'])

data.to_csv("C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Plot_files\\Original_france_annotation_exon_hist.tsv",
            sep='\t', index=True, header=True)