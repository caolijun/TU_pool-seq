#!/bin/bash
#SBATCH -J clj_pi
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 10
#SBATCH -o %j.out
#SBATCH -e %j.err

samples=pop1
poosize=400

prog=/work/share/kdy_caolijun/software/popoolation_1.2.2/Variance-sliding.pl
#mkdir piDscan
#mkdir piDscan/pi
# Measure pi for the populations
mpileup=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pileups/subsampled_depth50_${samples}.mpileup
out=/work/home/caolijun/tu_pool_seq/2popoolation/scripts/pi/${samples}_subsampled_depth50_pi_5k5k

perl $prog --input $mpileup --output $out --measure pi --window-size 5000 --step-size 5000 --min-count 1 --min-coverage 10 --min-qual 20 --pool-size $poosize --fastq-type sanger


# prepare data for plot
awk -F" " '{print $1,$2,$5}' $out > ${out}_1
sed -i 's/chr0//g' ${out}_1
awk '{ print FNR " " $0}' ${out}_1 > ${out}_ID
sed -i '1i\ID CHROM POS D'  ${out}_ID  ##加上表头,这里是空格:
sed -i 's/ /	/g' ${out}_ID ##将空格替换为tab
sed -ie '/na/d' ${out}_ID ##将na行删掉
rm ${out}_1 ${out}_IDe

