#!/bin/bash

# fastQC is for the quality control of raw data.

# First, be sure about have all the software from the previous requirements installed.
# Now, put all the fastq files with the raw data in the same directory.

# Create a new one to put in the results.
mkdir ./../results/FastQC_Reports

# Run FastQC for the analysis
for i in ~/../data/*.fastq;
do fastqc $i;
done

# Move your data to a different directory
mv *.zip *.html ./../results/FastQC_Reports/

# You need to examinate every one individually as is described on:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
