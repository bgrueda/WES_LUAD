# Workflow for the genomic variants identification

This workflow allows to identify genomic variants (germline and somatic) with some easy following steps. 

#### 1_Quality.sh

We use **FastQC** for this first quality analysis of the *fasta* files (raw data).

```bash
for i in *.fastq;
do fastqc $i;
done
```

once the report files were created, move them to a separate directory with: 

```bash
mkdir FastQC_Reports
mv *.zip *.html FastQC_Reports/
```

our report will contain some graphics with a colored signaling about:  

+ Basic statistics
+ Per base sequence quality
+ Per sequence quality scores
+ Per base sequence content
+ Per base GC content
+ Per sequence GC content 
+ Per base N content
+ Sequence Lenght Distribution
+ Sequence Duplication Levels
+ Overrepresented sequences 
+ Kmer Content

More details in: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

![Fastqc](https://github.com/bgrueda/WES_LUAD/blob/main/figures/Fastqc_ex.png)

This is an example of a good quality report (Per base sequence quality)

#### 2_Correction.sh

Once we have reviewed the detail of the quality of the reads, if necessary, we can make different types of corrections using **Trimmomatic**.

```bash
for i in *R1.fastq;
do java -jar /route/trimmomatic-0.39.jar PE -phred33 $i ${i%?.fastq}2.fastq trimmed_$i unpaired_$1 trimmed${i%?.fastq}2.fastq unpaired${i%?.fastq}2.fastq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:25 MINLEN:50;
done
```

The options should be chosen according to the state of the quality of the readings. It is highly recommended to check the changes made with a new fastqc analysis.  



#### 3_Mapping.sh

For the alignment to a reference genome we use **bwa** and a *fasta* file of the reference genome. We firstly needs to create an index or download it from the same web page source. 

The Broad Institute's genomics public data is available in: 

https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false 

```bash
bwa index hg38.fasta
for i in trimmed*R1.fastq;
do bwa mem hg38.fasta $i ${i%1.fastq}2.fastq > aln_${i%R1.fastq}.sam 2>bwa_${i%R1.fastq}.log;
done &
```

The output of this tool is a *SAM file*. 



#### 4_Preprocessing.sh

This step is divided in three parts:

+ Depth and coverage
+ Mark duplicates
+ BQSR

The first thing to do is changing the *SAM file* to its binary *BAM file* in order to have more efficient management. 

```bash
for i in aln*;
do samtools view -Sb $i > ${i%.sam}.bam 2>sb${i%.sam}.log | samtools sort ${i%.sam}.bam sort_${i%.sam}.bam 2>sort${i%.sam}.log;
done
```



*a. depth and coverage*

For this kind of analysis we use the **Samtools** suite.

```bash
for i in *.bam;
do samtools depth -a $i > ../results/depth/depth_${i%.bam}.txt 2>d_${i%.sam}.log;
samtools coverage $i -o ../results/coverage/cov_${i%.bam}.txt 2>c_${i%.sam}.log;
samtools coverage -m $i -o ../results/coverage/hist_${i%.bam}.txt 2>h_${i%.sam}.log;
done
```

The content of the files are: 

+ *depth.txt*

  ```bash
  ==> st_depth_nozero.txt <==
  chr1	825421	1
  chr1	825422	1
  chr1	825423	1
  chr1	825424	1
  chr1	825425	1
  chr1	825426	1
  chr1	825427	1
  chr1	825428	1
  chr1	825429	1
  chr1	825430	1
  ```

+ *cov.txt*

  ```bash
  ==> st_cov_data.txt <==
  #rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
  chr1	1	248956422	9654	753636	0.302718	0.00380098	66	58.4
  chr2	1	242193529	6476	508264	0.209859	0.00260629	66.3	58.9
  chr3	1	198295559	5694	448883	0.226371	0.00281306	66.3	59.4
  chr4	1	190214555	3575	280702	0.147571	0.00182409	66.4	58.5
  chr5	1	181538259	4068	319171	0.175815	0.00219193	66.3	58.8
  chr6	1	170805979	4112	322905	0.189048	0.00235587	66.3	56.5
  chr7	1	159345973	4492	356067	0.223455	0.0027579	66.3	57.3
  chr8	1	145138636	3241	256371	0.176639	0.00217759	66.2	57.2
  chr9	1	138394717	3881	306134	0.221204	0.00275086	65.9	58.6
  chr10	1	133797422	3892	306022	0.22872	0.002845	66.2	58.4
  
  ```

+ *hist.txt*

  ```bash
  ==> st_cov.txt <==
  chr1 (248.96Mbp)
  >   0.87% │                               █                  │ Number of reads: 9654
  >   0.77% │                              ██                  │     (1 filtered)
  >   0.68% │   ▅    ▁                     ██                  │ Covered bases:   753.6Kbp
  >   0.58% │   █ ▄  █                     ██                  │ Percent covered: 0.3027%
  >   0.48% │▄ ▁███  █                     ██▇       ▅▃        │ Mean coverage:   0.0038x
  >   0.39% │█▃████▆▇█▁            ▂       ███       ██   ▄    │ Mean baseQ:      66
  >   0.29% │██████████▇           █       ███       ██   █▂   │ Mean mapQ:       58.4
  >   0.19% │███████████       ▃  ▃█      ▆███▄▄▁▆   ██▂  ███ ▆│ 
  >   0.10% │███████████▇▇▆▁▄ ▄█ ▃███▁    ████████▁ ▁███▄▇███▂█│ Histo bin width: 4.98Mbp
  >   0.00% │████████████████▆██▅█████▁  ▁█████████▄███████████│ Histo max bin:   0.96766%
            1       49.79M    99.58M   149.37M   199.17M    248.96M 
  ```

  

For the depth files you can construct a plot to better visualization.

```R
setwd("~/bin")

# Call the file with the depth information from samtools depth command.
tab <- read.table ("depth.txt", 
                   header= F)
attach (tab)
names (tab)

# The plot will be saved as pdf
pdf("~/../results/depth1.pdf", width=6, height=3)

# Set up the columns of coordinate (x) and depth (y) values for the plot
x <- c(V2)
y <- c(V3)

# Make the plot
plot (x,y, 
      main="Depth",
      col= 'navy', 
      type= "l")

dev.off()
```

![Rplot](https://github.com/bgrueda/WES_LUAD/blob/main/figures/depth.jpeg)



*b. Mark duplicates*

The duplicates are ...To identify the duplicates we use the tool *MarkDuplicates* from **Picard**.

```bash
for i in sort*.bam;
do java -jar picard.jar MarkDuplicates I=$i O=md_$i M=metrics_$i.txt 2>md_${i%.bam}.log;
done
```



*c. BQSR (Base Quality Score Recalibration)*

For this recalibration 

```bash
for i in md_*;
do gatk BaseRecalibrator -I $i -R hg38.fasta --known-sites snp_hg38.vcf -O br_$i 2>md_${i%.txt}.log; 
done

for i in br_*;
do gatk ApplyBQSR -R hg38.fasta -I $i --bqsr-recal-file xBQSR_table_$i -O $i-BQSR.bam; 
done
```



#### 5_Somatic_Var.sh

```bash
for i in BQSR*; 
do gatk Mutect2 -R hg38.fasta -I ${i%BQSR_} -O ${i%BQSR_}.vcf 2>Mut${i%BQSR_}.log;
done
```



##### *Annotation*

For the visualization of the annotated vcf file, we run the next script. 

*maf_vis.R*

```bash
setwd("~/bin")

library(maftools)

# To call the data from the vcf file
xann = read.maf(maf = "file.maf")
xann

# To make the summary of the data in a separate table
getSampleSummary(xann)
getGeneSummary(xann)
getFields(xann)
write.mafSummary(maf = xann, basename = 'xann')

# To make a summary plot of the vcf data
plotmafSummary(maf = xann, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE + scale_color_viridis(discrete = TRUE, option = "D"))
```

```bash
==> geneSummary.txt
Hugo_Symbol Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins Missense_Mutation Nonsense_Mutation Nonstop_Mutation Splice_Site total MutatedSamples AlteredSamples
MUC16	0	0	0	0	9	0	0	0	9	1	1
TTN	0	1	0	0	5	0	0	1	7	1	1
KLC1	1	0	0	0	4	0	0	1	6	1	1

==> sampleSummary.txt
Tumor_Sampe_Barcode	Frame_Shift_Del	Frame_Shift_Ins	In_Frame_Del	In_Frame_Ins	Missense_Mutation	Nonsense_Mutation	Nonstop_Mutation	Splice_Site	total
__UNKNOWN__	20	63	15	25	1264	51	3	65	1506

==> summary.txt
ID	summary	Mean	Median
NCBI_Build	hg38	NA	NA
Center	__UNKNOWN__	NA	NA
Samples	1	NA	NA
nGenes	1236	NA	NA
Frame_Shift_Del	20	20	20
Frame_Shift_Ins	63	63	63
In_Frame_Del	15	15	15
In_Frame_Ins	25	25	25
Missense_Mutation	1264	1264	1264
Nonsense_Mutation	51	51	51
Nonstop_Mutation	3	3	3
Splice_Site	65	65	65
total	1506	1506	1506
```

![Image](https://github.com/bgrueda/WES_LUAD/blob/main/figures/summary_vcf.jpeg)





