# Whole Exome Sequencing of Lung Adenocarcinoma
This is a workflow to analyze WES (*Whole exome sequencing*) data in therms of the genomic variants identification.

## Motivation
Since cancer recently have been defined as a genome disease due its genomic instability, knowing a more complete landscape of mutations became one of the most important aspects to understand the particular carcinogenic process of each tumor in a more precise way. Likewise, identifying similarities and differences in relation to subtype, risk factors and even associated with populations is crucial. particular tumor, as well as differences between types, risk factors and even the origin of the population.
This could lead to better treatments using target drugs as an approach to the adoption of precision medicine.

## Workflow
This image below shows the steps of the preprocessing and processing of the data  for the identification of genomic variants and the software needed for this purpose.

![IMAGE](https://github.com/bgrueda/WES_LUAD/blob/main/Workflow.jpg)

## Content
This repository has the following organization:

```
archive

bin
  1_fastqc.sh
  2_trimmomatic.sh
  3_bwa.sh
  4_Preprocessing.sh

data
  data.md
  
figures

metadata
```

**archive(/archive)**. This directory have all the non-ready scripts.

**[bin](/bin)**. All the scripts needed for the analysis.

**[data](/data)**. This contains all the links for the data files that I used for this workflow.

**[figures](/figures)**. This directory contains all the figures that resulted from this analysis.

**metadata(/metadata)**. In this you can find all the tables, figures that are not results but are related to this analysis.



## Data specifications
**Data files**
1. The data that are needed for this kind of analysis are **fastq files** which are the result of sequencing of your samples (raw data).
2. A reference genome version (GRCh38). This can be found in one of the links in *data* folder.

**Sample preparation requirements**
+ **DNA** samples
+ *Tumor and/or normal* tissue
+ *Exon enrichment*
+ Sequenced in a NGS platform (Illumina Hiseq 2000, paired ends)
+ Calculating an average depth of **30X** approximately.  

## Requirements and software versions
- [FastQC (v.0.11.9)](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
- [Trimmomatic (v.0.39)](http://www.usadellab.org/cms/?page=trimmomatic)
- [BWA (v.0.7.12-r1039)](https://sourceforge.net/projects/bio-bwa/)
- [GATK (v.4.1.9.0)](https://github.com/broadinstitute/gatk/releases)
- [Samtools (v.1.10)](https://sourceforge.net/projects/samtools/)
- [R (v.3.4.2)](https://cran.r-project.org)
