#!/bin/bash

# Samtools (v.1.10) have some tools to manage the sam and bam files and do the depth and coverage analysis.
# GATK (v.4.1.9.0) and Picard (v.2.22.8) to mark duplicates and do the base recalibration.

# Convert sam file to bam format in order to have a more efficient management.
for i in aln*;
do samtools view -Sb $i > ${i%.sam}.bam 2>sb${i%.sam}.log | samtools sort ${i%.sam}.bam sort_${i%.sam}.bam 2>sort${i%.sam}.log;
done

# Depth & Coverage
for i in *.bam;
do samtools depth -a $i > ../results/depth/depth_${i%.bam}.txt 2>d_${i%.sam}.log;
samtools coverage $i -o ../results/coverage/cov_${i%.bam}.txt 2>c_${i%.sam}.log;
samtools coverage -m $i -o ../results/coverage/hist_${i%.bam}.txt 2>h_${i%.sam}.log;
done

# Mark duplicates in the bam files
for i in sort*.bam;
do java -jar picard.jar MarkDuplicates I=$i O=md_$i M=metrics_$i.txt 2>md_${i%.bam}.log;
done

#  Base Quality Score Recalibration
# You could use the vcf file from the Broad Institute or another that you know, even you could use more than one.
for i in md_*;
do gatk BaseRecalibrator -I $i -R hg38.fasta --known-sites snp_hg38.vcf -O br_$i 2>md_${i%.txt}.log;
done

# With this step, we are ready to continue to the variant calling step. Congratulations!
