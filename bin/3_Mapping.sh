#!/bin/bash

# BWA (v.0.7.12-r1039) align the reads to a reference genome
# You need to download the reference genome version that you prefer
# In this case,we are going to use: GRCh38
# It can be downloaded from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/

# Firstly we need to create an index
bwa index ./../data/hg38.fasta
mv *.fai ./../data/

# Do the alignment with the trimmed files.
for i in ./../results/2_Correction/trimmed*R1.fastq;
do bwa mem ./../data/hg38.fasta $i ${i%1.fastq}2.fastq > ./../results/3_Mapping/aln_${i%R1.fastq}.sam 2>bwa_${i%R1.fastq}.log;
done &
