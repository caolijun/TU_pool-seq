#!/bin/bash
#SBATCH -J clj_treemix
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -o %j.out
#SBATCH -e %j.err
module load compiler/gcc/10.2.0
input=/work/home/caolijun/tu_pool_seq/2popoolation/treemix/tu_22pops/treemix_input_22pops_subsample.gz
treemix=/work/share/kdy_caolijun/software/treemix-1.13/bin/treemix

# Running Treemix using SCCD as root; no migration

for j in {1..1000}
    do
   $treemix -i $input -bootstrap -k 100 -m 0 -root SCCD -o mite_tree_1mig_bootstrap_${j}
done