#!/bin/bash
input=$1
ref_address="/mnt/h/NAFLD/out2/hg38"
while read -r line 
do 
        #cutadapt -j 0 -a TGGAATTCTCGGGTGCCAAGG --length 30 -o output_2/${line}/${line}.trimmed.fastq.gz output_2/${line}/${line}.fastq.gz
		#fastqc -t 24 -o ./result/fastqc output_lnc/${line}/${line}.trimmed.fastq.gz
        #fastqc -t 48 -o ./result_lnc/fastqc output_lnc/${line}/${line}_*.fastq.gz
        #hisat2 -p 48 -x ${ref_address}/genome -1 output_lnc/${line}/${line}_1.fastq.gz -2 output_lnc/${line}/${line}_2.fastq.gz -S ./result_lnc/bam/${line}.sam
        #samtools sort -@ 10 -o ./result_lnc/bam/${line}.bam ./result_lnc/bam/${line}.sam
        #rm ./result_lnc/bam/${line}.sam
        #java -jar -XX:ParallelGCThreads=8 picard.jar AddOrReplaceReadGroups I=result_lnc/bam/${line}.bam O=result_lnc/bam/${line}.RG.bam RGID=rg1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${line}
        #java -jar -XX:ParallelGCThreads=8 picard.jar MarkDuplicates I=result_lnc/bam/${line}.RG.bam M=result_lnc/bam/${line}.metrics.txt ASO=queryname O=result_lnc/bam/${line}.RG.dedup.bam
        #java -jar -XX:ParallelGCThreads=8 picard.jar SortSam I=result_lnc/bam/${line}.RG.dedup.bam O=result_lnc/bam/${line}.RG.dedup.sorted.bam SORT_ORDER=coordinate
        #samtools index -@ 8 ./result_lnc/bam/${line}.RG.dedup.sorted.bam ./result_lnc/bam/${line}.RG.dedup.sorted.bai
        #java -jar -Xmx20g -Xms10g -XX:ParallelGCThreads=8 gatk-package-4.2.6.1-local.jar SplitNCigarReads -R hg38/Homo_sapiens_assembly38.fasta -I result/bam/${line}.RG.dedup.sorted.bam -O result/bam/${line}.RG.dedup.sorted.sc.bam
        #samtools index -@ 8 ./result/bam/${line}.RG.dedup.sorted.sc.bam ./result/bam/${line}.RG.dedup.sorted.sc.bai
		#stringtie ./result_lnc/bam/${line}.RG.dedup.sorted.bam -e -B -G ${ref_address}/gencode.v46.annotation.gtf -o ./result_lnc/gtf/${line}.gtf -p 48 -A ./result_lnc/count/${line}
		#python ../../prepDE.py3 -i . -l 150
		#conda activate /mnt/k/CIRIquant_env/
		CIRIquant -t 48 -1 output_lnc/${line}/${line}_1.fastq.gz -2 output_lnc/${line}/${line}_2.fastq.gz --config circquant.yml -o ./results_circRNA/${line} -p ${line}
done <  ${input}
