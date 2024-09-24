#!/bin/bash
#SBATCH -J clj_unmerge
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -o %j.out
#SBATCH -e %j.err

popoolation2=/work/share/kdy_caolijun/software/popoolation2_1201

### Calculate allele frequencies
perl ${popoolation2}/snp-frequency-diff.pl --input ${input}/unmerged_22pops.mpileup_depth50.sync --output-prefix ${output}/tu_22pops_subsample_depth50.freq --min-count 30 --min-coverage 4 --max-coverage 100000

### Extract allele frequencies (major,minor) from rc file using a custom Python script
python rc_to_pop_v3.py ${output}/tu_22pops_subsample_depth50.freq_rc > ${output}/tu_22pops_subsample_depth50.freq_majmin
#
## Attach header with population names, same order as they were merged.
echo -e "chr\tpos\tLabS\tNMHH\tBJWDY\tLabR\tZJHZ2\tHNCS\tBJDX2\tGXNN\tQHHD\tZJHZ1\tLNDD\tYNYX\tHNHK\tHNCS\tJXNC\tSCCD\tSHPD\tSDRZ\tAHHN\tSXYQ\tSDQZ\tSDSG3" | cat - ${output}/tu_22pops_subsample_depth50.freq_majmin > ${output}/tu_22pops_subsample_depth50.freq_majmin_labeled
#
## Remove the first two columns (chromosome; position), leaving only one column per population
cut -f 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 ${output}/tu_22pops_subsample_depth50.freq_majmin_labeled > ${output}/treemix_input_22pops_subsample
