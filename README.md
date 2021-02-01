# Whole Exome Sequencing of Lung Adenocarcinoma
This is a workflow to analyze WES (*Whole exome sequencing*) data in therms of the genomic variants identification.

## Motivation

Lung cancer is one of the most common types of cancer and its importance is associated to the high mortality of this disease. Lung adenocarcinoma is the most common subtype of lung cancer.

![IMAGE](https://github.com/bgrueda/WES_LUAD/blob/main/figures/lung_cancer.jpg)

Since cancer recently have been defined as a genome disease due its genomic instability, knowing a more complete landscape of mutations became one of the most important aspects to understand the particular carcinogenic process of each tumor in a more precise way.

Likewise, identifying similarities and differences in relation to subtype, risk factors and even associated with populations is crucial. particular tumor, as well as differences between types, risk factors and even the origin of the population.

This could lead to better treatments using target drugs as an approach to the adoption of precision medicine.

## Workflow
This image below shows the steps of the preprocessing and processing of the data  for the identification of genomic variants and the scripts needed for this purpose.

![IMAGE](https://github.com/bgrueda/WES_LUAD/blob/main/figures/Workflow.jpg)

## Content
This repository has the following organization:

```
|--- bin
|   |-- 1_Quality.sh
|   |-- 2_Correction.sh
|   |-- 3_Mapping.sh
|   |-- 4_Preprocessing.sh
|   |-- 4.1_depth.R
|   |-- 5_Somatic_Var.sh
|   |__ 5.1_maf_vis.R
|
|--- data
|   |-- data.md
|  
|--- figures
|   |-- figures.md
|
|___ results
    |-- 1_Quality
    |    |__ fastqc_report.html
    |
    |-- 2_Correction
    |    |__ correction.md
    |
    |-- 3_Mapping
    |    |-- aln_example.sam
    |    |__ sort_example.bam
    |    
    |-- 4_Preprocessing
    |    |-- BQSR
    |    |    |-- BQSR_example.bai
    |    |    |-- BQSR_example.bam
    |    |    |-- BQSR_example.table
    |    |    |-- BQSR_example.table_addline.txt
    |    |    |__ br_example.bam
    |    |    
    |    |-- coverage_and_depth
    |    |    |-- cov_example.txt
    |    |    |-- depth_example.txt
    |    |    |__ hist_example.txt
    |    |    
    |    |__ duplicates
    |         |--md_example_metrics.txt
    |         |__md_example.bam
    |      
    |__ 5_Somatic_Var
         |__ sample.maf
         |__ sample.vcf

```

**[Workflow_Tutorial.md](https://github.com/bgrueda/WES_LUAD/blob/main/Workflow_Tutorial.md)** This file contains all the information for the use of the scripts in this repository.  

**[bin](/bin)**. All the scripts needed for the analysis.

**[data](/data)**. This contains all the links for the data files that I used for this workflow.

**[figures](/figures)**. This directory contains all the images that are related or resulted from this analysis.

**[results](/results)**. This directory contains all the results from the pipeline.

## Data specifications
**Data files**
1. The data that are needed for this kind of analysis are **fastq files** which are the result of sequencing of your samples (raw data).
2. A reference genome version (hg38). This can be found in one of the links in *data* folder.

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
    + [MAFtoos (v.2.6.0)](https://bioconductor.org/packages/release/bioc/html/maftools.html)

## References
+ Genomic Landscape of Non-Small Cell Lung Cancer in Smokers and Never-Smokers
Ramaswamy Govindan, Li Ding, Malachi Griffith, Janakiraman Subramanian, Nathan D. Dees, Krishna L. Kanchi, Christopher A. Maher, Robert Fulton, Lucinda Fulton, John Wallis, … Richard K. Wilson
Cell (2012-09) https://doi.org/f39hs4
DOI: 10.1016/j.cell.2012.08.024 · PMID: 22980976 · PMCID: PMC3656590  

+ Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines
Kyle Ellrott, Matthew H. Bailey, Gordon Saksena, Kyle R. Covington, Cyriac Kandoth, Chip Stewart, Julian Hess, Singer Ma, Kami E. Chiotti, Michael McLellan, … Armaz Mariamidze
Cell Systems (2018-03) https://doi.org/gf9twn
DOI: 10.1016/j.cels.2018.03.002 · PMID: 29596782 · PMCID: PMC6075717

+ Some more information is available in: https://gatk.broadinstitute.org/hc/en-us
