
#!/bin/bash

# module load it/aligners/blast/2.10.1+
# module load it/tools/bedtools/2.27.1


chromosome=$1

ref_fasta=/projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/PN40024.v4_11_05_21/PN40024.v4.REF.fasta
my_ref_fasta=/projects/novabreed/share/mliva/hifiasm_assembly/pn40024/PN40024_v5.1/PN40024_v5.1.fasta
gff_to_read=/projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/PN40024.v4_11_05_21/REF_annotations/PN40024.v4.1.REF.gff3
gff_to_Unchr=/projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/gff_translation/PN40024.v4.1.REF.translated.PN40024_v5.1.unchr.gff3
stat=/projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/gff_translation/stats.txt


gff_to_write=/projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/gff_translation/PN40024.v4.1.translated.v5.1_${chromosome}_.gff3

previous_chr=None
mRNA_start=None
mRNA_end=None
mRNA_row_p1=None
mRNA_row_p2=None
gene_row_p1=None
gene_row_p2=None
gene_start=None
gene_end=None
cnt=0
num_genes_traslated=0
num_genes_seen=0

while read chr source type start end score strand phase attr; do

    # if row starts with ## it will be skiped
    if [[ ${chr} == \#\#* ]]
        then
        # echo "Row's beginning is strange skipping"
        continue

    elif [[ ${chr} != ${chromosome} ]]
        then
        continue
    

    # if the row has a gene or an mRMA type, informations are stored in a variable for later use
    elif [[ ${type} == "gene" ]]
        then
        if [[ ${gene_start} != "None" ]]
            then
            # echo "Here i am in a gene row and i am writing"
            echo -e "${gene_row_p1}\tgene\t${gene_start}\t${gene_end}${gene_row_p2}" >> ${gff_to_write}
            num_genes_traslated=$(echo -e "${num_genes_traslated}+1" | bc )
            gene_start=None
            gene_end=None

        fi
        num_genes_seen=$(echo -e "${num_genes_seen}+1" | bc )
        gene_row_p2=$(echo -e "\t${score}\t${strand}\t${phase}\t${attr}")
        gene_row_p1=$(echo -e "${chr}\t${source}")

    elif [[ ${type} == "mRNA" ]]
        then
        if [[ ${mRNA_start} != "None" ]] 
            then

            # echo "Here I am in mRNA row and i am writing"
            echo -e "${mRNA_row_p1}\tmRNA\t${mRNA_start}\t${mRNA_end}${mRNA_row_p2}" >>  ${gff_to_write}

            if [[ ${gene_start} == "None" ]]
                then
                echo "###" >> ${gff_to_write}
            fi

            # echo "mrna strat: ${mRNA_start}"
            # echo "mrna end: ${mRNA_end}"
            mRNA_start=None
            mRNA_end=None
        
        fi
        mRNA_row_p2=$(echo -e "\t${score}\t${strand}\t${phase}\t${attr}")
        mRNA_row_p1=$(echo -e "${chr}\t${source}")

    else

        bed_start=$(echo -e "${start}-1" | bc )

        echo -e "${chr}\t${bed_start}\t${end}" | bedtools getfasta -fi ${ref_fasta} -bed - > seq_${chromosome}.fasta

        if [[ ${chr} != ${previous_chr} ]]
            then

            # echo "Extracting my reference sequence"

            /projects/novabreed/share/gmagris/software/seqkit grep ${my_ref_fasta} -p "${chr}" -o filter_ref_${chr}.fasta
            previous_chr=${chr}
        fi 

        if [[ ${chr} == "chrUn" ]]
            then
            /iga/scripts/gmagris/run_blastn -q seq_${chromosome}.fasta -d ${my_ref_fasta} -o /projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/gff_translation/tmp_file${chromosome} -t 2

            awk -v start=${start} '{if ($3=="sacc") next } {if ( $13 >= 99.7 && $10 >= $2*0.9 ) print $3"\t"$8"\t"$9"\t"sqrt((start-$8)^2)}' tmp_file${chromosome}.tab | head -n 1 > out_blast${chromosome}.txt

            # echo "Unknown chr has been blasted"

        else

            # echo "Blasting..."
            /iga/scripts/gmagris/run_blastn -q seq_${chromosome}.fasta -d filter_ref_${chr}.fasta -o /projects/assembly_long_reads/vitis_vinifera/PN40024/gene_prediction/gff_translation/tmp_file${chromosome} -t 2


            awk -v start=${start} '{if ($3=="sacc") next } {if ( $13 >= 99.5 && $10 >= $2*0.88 ) print $3"\t"$8"\t"$9"\t"sqrt((start-$8)^2)}' tmp_file${chromosome}.tab | sort -V -k 4 | head -n 1 > out_blast${chromosome}.txt


        fi

        empy_file=$(wc -l out_blast${chromosome}.txt)
        # echo "${empy_file}"

        if [[ $empy_file == 0* ]]
        
            then 
            # echo "No rows in  the blast output file, skip row"
            #### Insert someting to count the number of bases that are not been found ###########
            cnt=$(echo -e "${cnt}+1" | bc )
            continue
        
        else

            num_chr=$(cut -f 1 out_blast${chromosome}.txt)
            file_start=$(cut -f 2 out_blast${chromosome}.txt)
            file_end=$(cut -f 3 out_blast${chromosome}.txt)

            if [[ ${gene_start} == "None" ]]
                then

                gene_start=${file_start}
                gene_end=${file_end}

            fi
            if [[ ${gene_start} != "None" ]]
                then
                if (( ${gene_start} > ${file_start} ))
                    then
                    gene_start=${file_start}

                fi

                if (( ${gene_end} < ${file_end} ))
                    then
                    gene_end=${file_end}

                fi
            fi

            if [[ ${mRNA_start} == "None" ]]
                then

                mRNA_end=${file_end}
                mRNA_start=${file_start}

            fi
            if [[ ${mRNA_start} != "None" ]]
                then
                if (( ${mRNA_start} > ${file_start} ))
                    then 
                    mRNA_start=${file_start}
                fi

                if (( ${mRNA_end} < ${file_end} ))
                    then
                    mRNA_end=${file_end}

                fi
            fi

            echo -e "${num_chr}\t${source}\t${type}\t${file_start}\t${file_end}\t${score}\t${strand}\t${phase}\t${attr}" >> ${gff_to_write}

        fi

    fi

done < ${gff_to_read}

echo -e "${chromosome}: ${cnt} lines skipped due too not finding any result during the blasting.\n: ${num_genes_seen} gene seen ${num_genes_traslated} gene traslocated\n" >> ${stat}
