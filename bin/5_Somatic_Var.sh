# Fot the somatic variants calling of one tumor sample.
# We need to have our reference genome.

mkdir ./../results/somatic_var

for i in ./../results/preprocessing/BQSR*;
do gatk Mutect2 -R hg38.fasta -I ${i%BQSR_} -O ./../results/somatic_var/${i%BQSR_}.vcf 2>Mut${i%BQSR_}.log;
done

# We can do an additional filtering of the obtained variants, however we need a vcf with the variants that we want to filter (contamination and/or others).
for i in ./../results/preprocessing/*.vcf;
do gatk FilterMutectCalls -R hg38.fasta -V $i --contamination-table contamination.table --tumor-segmentation segments.tsv -O ./../results/somatic_var/Filt_$i 2>Filt$i.log;
done
