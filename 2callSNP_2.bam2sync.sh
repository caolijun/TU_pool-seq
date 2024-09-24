#!/bin/bash
#SBATCH -J CLJ_pop
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 30
#SBATCH -o %j.out
#SBATCH -e %j.err

genome=/public/home/caolijun/tu_pool_seq/2popoolation/genome/tu.genome.ipm_v2.fasta
out=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pileups

find /work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/ -name *.sorted.rmdup.bam > bam.list

### Step1. run parallelly for each partition, then run subsequent steps
samtools mpileup -Q 20 -r chr01:1-17000000 -b bam.list -f $genome > $out/rawpile_all_chr01_1.mpileup
samtools mpileup -Q 20 -r chr01:17000001-33567794 -b bam.list -f $genome > $out/rawpile_all_chr01_2.mpileup
samtools mpileup -Q 20 -r chr02:1-15000000 -b bam.list -f $genome > $out/rawpile_all_chr02_1.mpileup
samtools mpileup -Q 20 -r chr02:15000001-29602728 -b bam.list -f $genome > $out/rawpile_all_chr02_2.mpileup
samtools mpileup -Q 20 -r chr03:1-12000000 -b bam.list -f $genome > $out/rawpile_all_chr03_1.mpileup
samtools mpileup -Q 20 -r chr03:12000001-24524128 -b bam.list -f $genome >  $out/rawpile_all_chr03_2.mpileup

# Step2. Remove indels
# Find indels
perl ${popoolation1_path}/basic-pipeline/identify-genomic-indel-regions.pl --input $out/rawpile_all.mpileup --output $out/indelregions_all.gtf --indel-window 5

# Remove indels
perl ${popoolation1_path}/basic-pipeline/filter-pileup-by-gtf.pl --gtf $out/indelregions_all.gtf --input $out/rawpile_all.mpileup --output $out/rmindels_all.mpileup

# Step3. Split mpileup
for i in {1..22}
do
col1=$(expr ${i} \* 3 + 1)
col2=$(expr ${i} \* 3 + 2)
col3=$(expr ${i} \* 3 + 3)
cut -d "	" -f1,2,3,${col1},${col2},${col3} rmindels_all.mpileup > rmindels_pop${i}.mpileup
#awk -F"\t" 'OFS="\t"{print $1,$2,$3,$'${col1}',$'${col2}',$'${col3}'}' rmindels_all.mpileup > rmindels_pop${i}.mpileup
done

# Step4. Subsample
for i in 10
do
perl ${popoolation1_path}/basic-pipeline/subsample-pileup.pl --input $out/rmindels_pop${i}.mpileup --output $out/subsampled_depth50_pop${i}.mpileup --method withoutreplace --fastq-type sanger --min-qual 20 --target-coverage 50 --max-coverage 5000
done

# Step5. mpileup2sync
for i in {1..22}
do
java -ea -Xmx7g -jar ${popoolation2}/mpileup2sync.jar --input $out/subsampled_depth50_pop${i}.mpileup --output $out/subsampled_depth50_pop${i}.mpileup.sync --fastq-type sanger --min-qual 20 --threads 30
done


