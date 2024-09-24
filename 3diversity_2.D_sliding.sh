#!/bin/bash
#SBATCH -J clj_D
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 10
#SBATCH -o %j.out
#SBATCH -e %j.err

samples=pop1
poosize=400
#samples2=SGLJC
#samples2_ps=400
mpileup=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pop2024/pileups/subsampled_depth50_pop1.mpileup
#mpileup=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pileups/subsampled_depth50_${samples}.mpileup
out=/work/home/caolijun/tu_pool_seq/2popoolation/scripts/TajimaD/${samples}_subsampled_depth50_pi_5k5k
prog=/work/share/kdy_caolijun/software/popoolation_1.2.2/Variance-sliding.pl
# Measure pi for the populations

perl $prog --input $mpileup --output $out --measure D --window-size 5000 --step-size 5000 --min-count 2 --min-coverage 10 --min-qual 20 --pool-size $poosize --fastq-type sanger

# prepare data for plot
awk -F" " '{print $1,$2,$5}' $out > ${out}_1
sed -i 's/chr0//g' ${out}_1
awk '{ print FNR " " $0}' ${out}_1 > ${out}_ID
sed -i '1i\ID CHROM POS D'  ${out}_ID  ##¼ÓÉÏ±íÍ·,ÕâÀïÊÇ¿Õ¸ñ:
sed -i 's/ /	/g' ${out}_ID ##½«¿Õ¸ñÌæ»»Îªtab
sed -ie '/na/d' ${out}_ID ##½«naÐÐÉ¾µô
rm ${out}_1 ${out}_IDe

