#!/bin/bash

# Trimmomatic is the tool that allows you to do a quality filtering of the raw data.

# This part coul change depending of the results of the fastQC
# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 -> Remove adapters
# LEADING:3 -> Remove leading low quality or N bases (below quality 3)
# TRAILING:3 -> Remove trailing low quality or N bases (below quality 3)
# SLIDINGWINDOW:4:15 -> Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
# MINLEN:50 -> Drop reads below the 50 bases long.
# CROP: 150 -> The number of bases to keep, from the start of the read.
# HEADCROP:15 -> The number of bases to remove from the start of the read.


for i in *R1.fastq;
do java -jar /route/trimmomatic-0.39.jar PE -phred33 $i ${i%?.fastq}2.fastq trimmed_$i unpaired_$1 trimmed${i%?.fastq}2.fastq unpaired${i%?.fastq}2.fastq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:25 MINLEN:50;
done

# Then is needed to run Fastqc again in order to check the quality of the files after the cleanning.

for i in trimmed*;
do fastqc $i;
done

mv *.zip *.html FastQC_Reports/
