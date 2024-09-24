#!/bin/bash
#SBATCH -J clj_unmerge
#SBATCH -p xhacnormalb
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -o %j.out
#SBATCH -e %j.err

samples1=pop1
samples2=pop2
samples3=pop3
samples4=pop4
samples5=pop5
samples6=pop6
samples7=pop7
samples8=pop8
samples9=pop9
samples10=pop10
samples11=pop11
samples12=pop12
samples13=pop13
samples14=pop14
samples15=pop15
samples16=pop16
samples17=pop17
samples18=pop18
samples19=pop19
samples20=pop20
samples21=pop21
samples22=pop22

popoolation2=/work/share/kdy_caolijun/software/popoolation2_1201
input=/work/home/caolijun/tu_pool_seq/2popoolation/tmp_pipe_data/pop2024/pileups
output=/work/home/caolijun/tu_pool_seq/2popoolation/treemix

## Filter 2 by 1
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/subsampled_depth50_${samples1}.mpileup.sync ${input}/subsampled_depth50_${samples2}.mpileup.sync > ${input}/unmerge_temp1.2
## Filter 3 by 1;2
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2 ${input}/subsampled_depth50_${samples3}.mpileup.sync > ${input}/unmerge_temp1.2.3
## Filter 4 by 1;2;3
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3 ${input}/subsampled_depth50_${samples4}.mpileup.sync > ${input}/unmerge_temp1.2.3.4
## Filter 5 by 1;2;3;4
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4 ${input}/subsampled_depth50_${samples5}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5
## Filter 6 by 1;2;3;4;5
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5 ${input}/subsampled_depth50_${samples6}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6
## Filter 7 by 1;2;3;4;5
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6 ${input}/subsampled_depth50_${samples7}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7
## Filter 8 by 1;2;3;4;5;6
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7 ${input}/subsampled_depth50_${samples8}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8
## Filter 9 by 1;2;3;4;5;6;7
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8 ${input}/subsampled_depth50_${samples9}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9
## Filter 10by 1;2;3;4;5;6;7;8
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9 ${input}/subsampled_depth50_${samples10}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10
## Filter 11 by 1;2;3;4;5;6;7;8;9
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10 ${input}/subsampled_depth50_${samples11}.mpileup.sync >  ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11
# Filter 12 by 1;2;3;4;5;6;7;8;9;10
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11 ${input}/subsampled_depth50_${samples12}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12 ${input}/subsampled_depth50_${samples13}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13 ${input}/subsampled_depth50_${samples14}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14 ${input}/subsampled_depth50_${samples15}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15 ${input}/subsampled_depth50_${samples16}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16 ${input}/subsampled_depth50_${samples17}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17 ${input}/subsampled_depth50_${samples18}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18 ${input}/subsampled_depth50_${samples19}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19 ${input}/subsampled_depth50_${samples20}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20 ${input}/subsampled_depth50_${samples21}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$4}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21 ${input}/subsampled_depth50_${samples22}.mpileup.sync > ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22
#### Throw away all remaining incomplete rows (keep inner join of dataset)
awk 'BEGIN {FS=OFS="\t"} {if(NF>4){print $0}}' ${input}/unmerge_temp1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22 > ${input}/unmerged_22pops_subsample_depth50.sync
#rm unmerge_temp*


