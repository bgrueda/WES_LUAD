#!/bin/bash

# Samtools (v.1.10) have some tools to manage the sam and bam files and do the depth and coverage analysis.
# GATK (v.4.1.9.0) and Picard (v.2.22.8) to mark duplicates and do the base recalibration.

# Convert sam file to bam format in order to have a more efficient management.
for i in ./../results/mapped/aln*;
do samtools view -Sb $i > ./../results/mapped/${i%.sam}.bam 2>sb${i%.sam}.log | samtools sort ${i%.sam}.bam ./../results/mapped/sort_${i%.sam}.bam 2>sort${i%.sam}.log;
done

# Depth & Coverage
mkdir ./../results/preprocessing

for i in ./../results/*.bam;
do samtools depth -a $i > ./../results/preprocessing/depth_${i%.bam}.txt 2>d_${i%.sam}.log;
samtools coverage $i -o ./../results/preprocessing/cov_${i%.bam}.txt 2>c_${i%.sam}.log;
samtools coverage -m $i -o ./../results/preprocessing/hist_${i%.bam}.txt 2>h_${i%.sam}.log;
done

# For the plot of the depth, run 4.1_depth.R script.

# Mark duplicates in the bam files
for i in ./../results/preprocessing/sort*.bam;
do java -jar picard.jar MarkDuplicates I=$i O=./../results/preprocessing/md_$i M=./../results/preprocessing/metrics_$i.txt 2>md_${i%.bam}.log;
done

#  Base Quality Score Recalibration
# You could use the vcf file from the Broad Institute or another that you know, even you could use more than one.
for i in ./../results/preprocessing/md_*;
do gatk BaseRecalibrator -I $i -R ./../data/hg38.fasta --known-sites ./../data/snp_hg38.vcf -O ./../results/preprocessing/br_$i 2>md_${i%.txt}.log;
done

# With this step, we are ready to continue to the variant calling step. Congratulations!
