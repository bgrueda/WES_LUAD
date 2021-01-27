## Looking for genomic variants

Since the NGS (Next Generation Sequencing) become more accesible, the possibilities opened up before us to know better the genomic landscape involved in some complex pathological conditions as cancer. 

This repository aims to bring together the different tools used on the identification of genomic variants at present, to generate a reproducible and easy-to-follow workflow.  

We present a pipeline based on **GATK** best practices with some other complementary tools to identificate somatic and germline (WIP) vriants from WES (whole exome sequencing) sample files. 



#### About the results...

This pipeline was probed with a single small sample (part of a complete lung adenocarcinoma sample sequenced). So the results obtained lack of a biological significant meaning. 

In the current state of this pipeline, it is possible to arrive at the annotation of somatic genomic variants. 

![](https://github.com/bgrueda/WES_LUAD/blob/main/figures/summary_vcf.jpeg)

Nonetheless, when we confirm a set of the resulting variants with the IGV from the Broad Institute, we found that the tested ones already exists.  

![](https://github.com/bgrueda/WES_LUAD/blob/main/figures/IGV.png)

Further analysis can be developed to get more meaningful interpretations and other clinical approaches (WIP).