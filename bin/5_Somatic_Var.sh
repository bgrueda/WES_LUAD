# Fot the somatic variants calling of one tumor sample.
# We need to have our reference genome.

mkdir ./../results/somatic_var

for i in ./../results/4_Preprocessing/BQSR*;
do gatk Mutect2 -R ./../data/hg38.fasta -I $i -O ./../results/5_Somatic_Var/${i%BQSR_}.vcf 2>Mut${i%BQSR_}.log;
done

# We can do an additional filtering of the obtained variants, however we need a vcf with the variants that we want to filter (contamination and/or others).
for i in ./../results/4_Preprocessing/*.vcf;
do gatk FilterMutectCalls -R ./../data/hg38.fasta -V $i --contamination-table ./../data/contamination.table --tumor-segmentation ./../data/segments.tsv -O ./../results/somatic_var/Filt_$i 2>Filt$i.log;
done

# For the annotation, we use Funcotator from GATK toolkit
# Firstly is needed to download the sources in data directory

gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download

for i in ./../results/4_Preprocessing/Filt_*;
do gatk Funcotator -R ./../data/hg38.fasta -V $i -O ./../results/5_Somatic_Var/$i.maf --output-file-format MAF --data-sources-path ./../data/SSources/ --ref-version hg38 2>Func$i.log;
done
