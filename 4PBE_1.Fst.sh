#!/bin/bash
#SBATCH -J CLJ_pop
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 30
#SBATCH -o %j.out
#SBATCH -e %j.err

popoolation1_path=/work/share/kdy_caolijun/software/popoolation_1.2.2
popoolation2=/work/share/kdy_caolijun/software/popoolation2_1201
sync=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pop2024/pileups/unmerged_22pops.mpileup_depth50.sync

perl parallel_fstsliding.pl 30 $sync tu_22pops 3000 400:600:400:400:800:600:600:200:100:800:500:300:600:600:600:400:600:600:100:600:400:500 5000 5000
