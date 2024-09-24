samples1=LabS
samples2=NMHH
samples3=BJDX6

#####get FST for each population pair
awk -F" " '{print $1,$2,$6}' tu_16pops_all.fst > ${samples1}_${samples2}_5k_ 50depth_3L.fst
sed -i 's/1:2=//g' ${samples1}_${samples2}_5k_ 50depth_3L.fst ##去掉1:2=

awk -F" " '{print $1,$2,$15}' tu_16pops_all.fst > ${samples1}_${samples3}_5k_ 50depth_3L.fst
sed -i 's/1:11=//g' ${samples1}_${samples3}_5k_ 50depth_3L.fst ##去掉1:2=

awk -F" " '{print $1,$2,$29}' tu_16pops_all.fst > ${samples2}_${samples3}_5k_ 50depth_3L.fst
sed -i 's/2:11=//g' ${samples2}_${samples3}_5k_ 50depth_3L.fst ##去掉1:2=

## Filter 2 by 1
awk 'BEGIN {FS=OFS=" "} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$3}' ${samples1}_${samples2}_5k_ 50depth_3L.fst ${samples1}_${samples3}_5k_ 50depth_3L.fst > temp1.2

## Filter 3 by 1;2
awk 'BEGIN {FS=OFS=" "} NR==FNR {a[$1,$2]=$0;next} {print a[$1,$2],$3}' temp1.2 ${samples2}_${samples3}_5k_ 50depth_3L.fst > temp1.2.3

### Throw away all remaining incomplete rows (keep inner join of dataset)
awk 'BEGIN {FS=OFS=" "} {if(NF>4){print $0}}' temp1.2.3 > ${samples1}_${samples2}_${samples3}_5k.fst


#head of ${samples1}_${samples2}_${samples3}.fst
#sed -i '1d' ${samples1}_${samples2}_${samples3}_5k.fst
sed -i '1i\CHROM POS S1S2 RS1 RS2' ${samples1}_${samples2}_${samples3}_5k.fst  ##加上表头
sed -i 's/ /	/g' ${samples1}_${samples2}_${samples3}_5k.fst 
sed -i '/na/d' ${samples1}_${samples2}_${samples3}_5k.fst
######calculate PBE and plot figs
Rscript manhattan_ggplot_5k_v1.R ${samples1}_${samples2}_${samples3}