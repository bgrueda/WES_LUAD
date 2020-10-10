# WES_LUAD
This is a workflow to analyze WES data in therms of genomic variants identification.

## Motivation 
Since cancer recently have been defined as a genome disease due its genomic instability, knowing a more complete lanscape of mutations became one of the most important aspects to undesrstand the particular carcinogenic process of the tumor.
This could lead to better treatments using target drugs as an approach to the adoption of precision medicine.

## Workflow 
This image below shows the steps of the preprocessing and processing of the datafor the identification of genomic variants and the software needed to this purpose.

[IMAGE] 

## Content
This repository has the following organization

## Data 
The data that are needed for this kind of analysis are **fastq files**.

**DNA** samples from *tumor and/or normal* tissue with an *exon enrichment* and sequenced in a NGS plataform, *paired ends* preferently. Calculating an average depth of **30X** approximadetly.  

## Requirements and software versions
- FastQC (v,0.11.9)
- Trimmomatic (v.0.39)
- BWA (0.7.12-r1039)
- GATK (v.4.1.8.0)
- Samtools (v.1.10)



