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

*depth and coverage*

```bash
for i in *.bam;
do samtools depth -a $i > ../results/depth/depth_${i%.bam}.txt 2>d_${i%.sam}.log;
samtools coverage $i -o ../results/coverage/cov_${i%.bam}.txt 2>c_${i%.sam}.log;
samtools coverage -m $i -o ../results/coverage/hist_${i%.bam}.txt 2>h_${i%.sam}.log;
done
```

The content of the files are: 

+ *depth.txt*

  ![depth](https://github.com/bgrueda/WES_LUAD/blob/main/figures/depth.png)

+ *cov.txt*

  ![cov](https://github.com/bgrueda/WES_LUAD/blob/main/figures/cov.png)

+ *hist.txt*

  ![depth](https://github.com/bgrueda/WES_LUAD/blob/main/figures/hist.png)



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

