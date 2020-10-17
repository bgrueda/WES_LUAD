#!/bin/bash

# Samtools (v.1.10) have some tools to manage the sam and bam files.

# Convert sam file to bam format in order to have a more efficient management.
for i in aln*;
do samtools view -Sb $i > ${i%.sam}.bam 2>sb${i%.sam}.log | samtools sort ${i%.sam}.bam sort_${i%.sam}.bam 2>sort${i%.sam}.log;
done

# GATK (v.4.1.8.0) and Picard (v.2.22.8) to mark duplicates
# MArk duplicates in the bam files
for i in sort*.bam;
do java -jar picard.jar MarkDuplicates I=$i O=md_$i M=metrics_$i.txt 2>md_$i.log;
done 
