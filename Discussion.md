# Looking for genomic variants

Since the NGS (Next Generation Sequencing) become more accesible, the possibilities opened up before us to know better the genomic landscape involved in some complex pathological conditions as cancer.

This repository aims to bring together the different tools used on the identification of genomic variants at present, to generate a reproducible and easy-to-follow workflow.  

We present a pipeline based on **GATK** best practices with some other complementary tools to identify *somatic* and *germline* (WIP) variants from WES (whole exome sequencing) sample files.  



### About the results...

This pipeline was probed with a single small sample (part of a complete lung adenocarcinoma sample sequenced), then I ran it on our server with the entire sample, corroborating the proposed workflow.

In the current state of this pipeline, it is possible to arrive at the annotation of somatic genomic variants.

![](https://github.com/bgrueda/WES_LUAD/blob/main/figures/summary_vcf.jpeg)

Fig1. Fragmented sample

![](https://github.com/bgrueda/WES_LUAD/blob/main/figures/summary_server.jpeg)

  Fig2. Complete sample.

  Further analysis can be developed to get more meaningful interpretations and other clinical approaches (WIP).
