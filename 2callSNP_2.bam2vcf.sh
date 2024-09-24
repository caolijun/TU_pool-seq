#!/bin/bash
#SBATCH -J bam2vcf
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 30
#SBATCH -o %j.out
#SBATCH -e %j.err

genome=/work/home/caolijun/tu_pool_seq/2popoolation/genome/tu.genome.ipm_v2.fasta
bamlist=bam.list

bcftools mpileup --threads 30 -f $genome -Ou --bam-list $bamlist | bcftools call --threads 30 -mv -Oz > tu_pool2_vcf.gz
