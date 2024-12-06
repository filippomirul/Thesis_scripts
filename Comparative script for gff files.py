# Comparative script for gff files
import csv
import datetime

"""
Comparatiove script to compare the PN reference annotation, coming from the translation, to the Funannotate annotation of PN.
"""

####

def tmp_file_writer(infile, outfile, gene_at_start):

    with open(infile, "r") as in_file, \
        open(outfile, "w", newline = "") as out_file:

        reader = csv.reader(in_file, dialect = "excel", delimiter = "\t")
        writer = csv.writer(out_file, dialect="excel", delimiter="\t")

        cc = 0

        # new file rows: gene\t ID\t start\t end\t strand\t num_transcripts\t num_non_overlapping_exons

        if gene_at_start:
            cnt = -1
            # rows: list containing list. Each list is representative of one gene and contains the information of each gene
            rows = []
            # exons_positions: list of list containing for each gene starting and ending of unique exons
            exons_positions = []
            chr_list = []
            chr_cnt = -1

            for row in reader:

                if len(row) != 9:
                    continue

                # if row[0] == "chrUn":
                #     continue

                # if row[0] == "SKIP":
                #     break

                if row[0] not in chr_list:
                    chr_cnt += 1
                    chr_list.append(row[0])
                    exons_positions.append([])

                if row[2] == "gene":
                    rows.append([])
                    cnt += 1
                    rows[cnt].append("gene")
                    rows[cnt].append(row[8].split(";")[0])  # ID
                    rows[cnt].append(row[3])    # start
                    rows[cnt].append(row[4])    # end
                    rows[cnt].append(row[6])    # stand
                    rows[cnt].append(0)   # num_transcripts
                    rows[cnt].append(0)   # num non overlapping exons

                if row[2] == "mRNA":
                    rows[cnt][5] += 1

                if row[2] == "exon":
                    if row[3] in exons_positions[chr_cnt] and row[4] in exons_positions[chr_cnt]:
                        continue
                    else:
                        rows[cnt][6] += 1
                        exons_positions[chr_cnt].append(row[3])
                        exons_positions[chr_cnt].append(row[4])
        else:

            cnt = 0
            rows = [[]]
            exons_positions = []
            chr_list = []
            chr_cnt = -1

            rows[cnt].append(None)  # gene
            rows[cnt].append(None)  # ID
            rows[cnt].append(None)  # start
            rows[cnt].append(None)  # end
            rows[cnt].append(None)  # strand
            rows[cnt].append(0)   # num_transcripts
            rows[cnt].append(0)   # num non overlapping exons

            for row in reader:

                # if row[0] == "SKIP":
                #     break

                # if row[0] =="chrUn":
                #     cc += 1
                #     continue

                if len(row[0]) != 5 or len(row) != 9 or len(row[3]) < 5 or len(row[4]) < 5:
                    cc += 1
                    continue

                if row[0] not in chr_list:
                    #print(row[0])
                    chr_list.append(row[0])
                    exons_positions.append([])
                    chr_cnt += 1

                if row[2] == "gene":
                    rows[cnt][0] = "gene"
                    rows[cnt][1] = row[8].split(";")[0] # ID
                    rows[cnt][2] = row[3]   # start
                    rows[cnt][3] = row[4]   # end
                    rows[cnt][4] = row[6]   # strand
                    rows.append([None, None, None, None, None, 0, 0])
                    cnt += 1

                if row[2] == "mRNA":
                    rows[cnt][5] += 1

                if row[2] == "exon":
                    if row[3] in exons_positions[chr_cnt] and row[4] in exons_positions[chr_cnt]:
                        cc += 1
                        continue

                    else:
                        rows[cnt][6] += 1
                        exons_positions[chr_cnt].append(row[3])
                        exons_positions[chr_cnt].append(row[4])

                cc += 1
                # print(f"Linea: {cc}")
                # print(f"Pre: {row[3], row[4]}")
                # print(f"Post: {int(row[3]), int(row[4])}")
        
        for line in rows:
            writer.writerow(line)
    
    return exons_positions

## Inferring stats from intermedite file

def stats_from_annotation_comparison(reference_tmp, alternative_tmp):

    with open(reference_tmp, "r") as reff_file, \
        open(alternative_tmp, "r") as my_ann_file:

        reader_1 = csv.reader(reff_file, dialect="excel", delimiter="\t")
        reader_2 = csv.reader(my_ann_file, dialect = "excel", delimiter = "\t")

        cnt_gene_reader_1 = 0
        cnt_gene_reader_2 = 0

        cnt_transcripts_reader_1 = 0
        cnt_transcripts_reader_2 = 0

        cnt_exons_reader_1 = 0
        cnt_exons_reader_2 = 0

        for row in reader_1:

            cnt_gene_reader_1 += 1
            cnt_exons_reader_1 += int(row[6])
            cnt_transcripts_reader_1 += int(row[5])

        for row in reader_2:

            cnt_gene_reader_2 += 1
            cnt_exons_reader_2 += int(row[6])
            cnt_transcripts_reader_2 += int(row[5])

    print(f"Reference annotation:\n")
    print(f"Number of gene: {cnt_gene_reader_1}")
    print(f"Number of total transcripts: {cnt_transcripts_reader_1}")
    print(f"Number of total exons: {cnt_exons_reader_1}")
    print(f"Number of transcripts per gene: {cnt_transcripts_reader_1/cnt_gene_reader_1}")
    print(f"Number of exon per gene: {cnt_exons_reader_1/cnt_gene_reader_1}\n")
    print(f"Alternative annotation:\n")
    print(f"Number of gene: {cnt_gene_reader_2}")
    print(f"Number of total transcripts: {cnt_transcripts_reader_2}")
    print(f"Number of total exons: {cnt_exons_reader_2}")
    print(f"Number of transcripts per gene: {cnt_transcripts_reader_2/cnt_gene_reader_2}")
    print(f"Number of exon per gene: {cnt_exons_reader_2/cnt_gene_reader_2}\n")

#### Comparing alternative annotation with reference

def check_nexts(number, exon_list, var = 30):
    
    numb = int(number)
    tt = 0

    for i in exon_list:
        if numb > (int(i) + var) or numb < (int(i) - var):
            continue
        else:
            tt = 1
        
    if tt == 1:
        return True
    else:
        return False

def comparing_exons(alternative, ref_exons_starts_ends):
    with open(alternative, "r") as file:

        file__ = csv.reader(file, dialect="excel", delimiter="\t")

        cnt_consensus = 0
        cnt_missensus = 0
        cnt_partial = 0

        chr_list = []
        chr_cnt = -1

        cc = 0

        for row in file__:

            if row[0] == "SKIP":
                break

            if row[0] == "chrUn":
                continue

            if len(row) != 9 or len(row[0]) != 5 or len(row[3]) <= 1 or len(row[4]) <= 1:
                continue

            if row[0] not in chr_list:
                chr_list.append(row[0])
                chr_cnt += 1

            if row[2] == "exon":

                if check_nexts(row[3], ref_exons_starts_ends[chr_cnt]) and check_nexts(row[4], ref_exons_starts_ends[chr_cnt]):
                    cnt_consensus += 1

                elif check_nexts(row[3], ref_exons_starts_ends[chr_cnt]) or check_nexts(row[4], ref_exons_starts_ends[chr_cnt]):
                    cnt_partial += 1

                else:
                    cnt_missensus += 1

                
                # if row[3] not in ref_exons_starts_ends[chr_cnt] and row[4] not in ref_exons_starts_ends[chr_cnt]:
                #     cnt_missensus += 1

                # if row[3] in ref_exons_starts_ends[chr_cnt] and row[4] in ref_exons_starts_ends[chr_cnt]:
                #     cnt_consensus += 1

                # else:
                #     cnt_partial += 1

            # print(int(row[3]), int(row[4]))
            # print(cc)
            cc += 1

    return (cnt_consensus, cnt_partial, cnt_missensus)

### Reference-translation comparison

def comparison_france_translated_REF(original_ref, traslated_ref):

    tot_number_exon_france_ref = 0

    with open(original_ref, "r") as file:

        file__ = csv.reader(file, dialect="excel", delimiter="\t")
        france_ref_exons = []
        tot_exons = 0

        for row in file__:

            if row[0] == "SKIP":
                break

            if row[0] == "chrUn":
                continue

            if len(row[0]) != 5 or len(row) != 9 or len(row[3]) < 5 or len(row[4]) < 5:
                continue

            if row[2] == "exon":
                tot_exons += 1
                france_ref_exons.append(row)
                tot_number_exon_france_ref += 1

    with open(traslated_ref, "r") as file:

        file__ = csv.reader(file, dialect="excel", delimiter="\t")
        traslated_exons = []

        for row in file__:

            # if row[0] == "SKIP":
            #     break

            # if row[0] == "chrUn":
            #     continue 

            if len(row[0]) != 5 or len(row) != 9 or len(row[3]) < 5 or len(row[4]) < 5:
                continue

            if row[2] == "exon":
                traslated_exons.append(row)

    cnt = 0
    tot_len = 0
    len_exon_un = []

    for i in france_ref_exons:
        if i not in traslated_exons:
            cnt += 1
            tot_len += int(i[4]) - int(i[3])
            len_exon_un.append(int(i[4]) - int(i[3]))

    print("Untraslated exons\n")
    print(f"Number of exons untraslated: {cnt}")
    print(f"Percentage of untraslated exons: {cnt/tot_exons}")
    print(f"Mean length of the untraslated exons: {tot_len/cnt}\n")

    return len_exon_un

########    ######

debug = False

## First part creates the intermediate files

reference_annotation_file = "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\merged_transl_PN_france.gff3"
reference_tmp_annotation_file = "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\tmp_merged_traslated_PN_france.gff"

alternative_annotation_file = "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\PN_noupdate_genome.gff3"
alternative_tmp_annotation_file = "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Annotations\\PN_files\\tmp_PN_noupdate_genome.gff3"


# origina_REF_france_gff = "C:\\Users\\filoa\\Desktop\\Programming_trials\\IGA_scripts\\Original_untraslated_REF_france.gff3"

reference_tuple = (reference_annotation_file, reference_tmp_annotation_file)
alternative_tuple = (alternative_annotation_file, alternative_tmp_annotation_file)

# Some gff files put the gene at the end some at the beginning
# this information is necessary to the correct extraction of information

# before run the comparison after the reference, so ref exon information are kept

print(f"[{datetime.datetime.now()}]")
print(f"File reference: {reference_annotation_file}")
print(f"File alternative annotation: {alternative_annotation_file}")

at_start = True
alternativi_exons_positions = tmp_file_writer(alternative_tuple[0], alternative_tuple[1], gene_at_start= at_start)

at_start = False
reference_exons_positions = tmp_file_writer(reference_tuple[0], reference_tuple[1], gene_at_start= at_start)

print(f"[{datetime.datetime.now()}]")
stats_from_annotation_comparison(reference_tmp_annotation_file, alternative_tmp_annotation_file)

if debug == True:
    print(alternativi_exons_positions[0][:20])
    print(reference_exons_positions[0][:40])

    print(len(alternativi_exons_positions))
    print(len(reference_exons_positions))

print(f"[{datetime.datetime.now()}]")
data = comparing_exons(alternative_annotation_file, reference_exons_positions)
# data = comparing_exons(reference_annotation_file, alternativi_exons_positions)

print("Alternative c'è in reference")
print(f"Consensus exons: {data[0]}")
print(f"Partial consensus exons: {data[1]}")
print(f"Non overlapping exons (FP): {data[2]}\n")
aa = data[0] + data[1]
bb = data[0] + data[1] + data[2]
spec = (data[0] + data[1]) / (data[0] + data[1] + data[2])
print(f"Specificity: {spec}\n")   # FP quelli del predict che non ci sono sul reference

print(f"[{datetime.datetime.now()}]")
data = comparing_exons(reference_annotation_file, alternativi_exons_positions)

print("Reference c'è in alternative")
print(f"Consensus exons: {data[0]}")
print(f"Partial consensus exons: {data[1]}")
print(f"Non overlapping exons (FN): {data[2]}\n")
sens = (data[0] + data[1]) / (data[0] + data[1] + data[2])
cc = data[0]+data[1]
dd = data[0] + data[1] + data[2]
print(f"Sensitivity: {sens}")   # FN reference che non ci sono nel predict
print(f"Accuracy: {(sens + spec) / 2}")
print(f"MY accuracy: {(aa+cc)/(bb+dd) }")

# print(f"[{datetime.datetime.now()}]")
# comparison_france_translated_REF(origina_REF_france_gff, reference_annotation_file)

################
#   OUTPUT
################

# [2024-09-24 17:34:33.385872]
# Reference annotation:

# Number of gene: 34596
# Number of total transcripts: 40719
# Number of total exons: 153390
# Number of transcripts per gene: 1.1769857787027402
# Number of exon per gene: 4.433749566423864

# Alternative annotation:

# Number of gene: 35734
# Number of total transcripts: 35734
# Number of total exons: 152429
# Number of transcripts per gene: 1.0
# Number of exon per gene: 4.265657357138859

# [2024-09-24 17:34:33.477569]
# Alternative c'è in reference
# Consensus exons: 95809
# Partial consensus exons: 22444
# Non overlapping exons (FP): 34176

# Specificity: 0.7757906959961687

# [2024-09-24 17:53:40.249915]
# Reference c'è in alternative
# Consensus exons: 125753
# Partial consensus exons: 43318
# Non overlapping exons (FN): 35908

# Sensitivity: 0.8248210792325067
# Accuracy: 0.8003058876143376
# MY accuracy: 0.8039103769361626