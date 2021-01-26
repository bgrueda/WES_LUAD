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