#!/bin/bash
# The shell script must be made executable by changing its permission bits. An example:
# 1 input values are the number of threads to use
# 2 input sync file 
# 3 output filename
# 4 max read depth
# 5 pool size
# 6 window size
# 7 step size
#$ chmod u+x *.sh
mkdir temp_popool_fstsliding
cd temp_popool_fstsliding
split -n l/$1 $2
# index the files
ls x* > files
# files2read=$(ls x*)
outfile=$3

# clean up the splits so that each one had a complete contig set within it
# files is the ordered list of these, so this can be used to walk through

count=1
# warning, this resets all numeric variables to being assigned to items of $files
array=($(ls x*))
file_count=$(wc -l files | cut -f1 -d ' ')
# echo $file_count
while read p; do
# 		if [ "$count" == 2 ]; then
# 			continue
# 		fi
		if [ "$count" = "$file_count" ]; then
			cat "$p" > $p.tail_blunted
			break
		fi
# 		echo $count
# 		echo $p
# 		echo ${array[$count]}
# 		following_file=${array[$count]}
# 		echo $following_file
		tail_leadingcontig=$(tail -1 $p | cut -f1) 
		head_followingcontig=$(head -1 ${array[$count]} | cut -f1)
# 		echo $tail_leadingcontig
# 		echo $head_followingcontig
		if [ "$tail_leadingcontig" = "$head_followingcontig" ]; then
			lines_head=$(grep "$head_followingcontig" -c ${array[$count]})
# 			echo $lines_head
			head -"$lines_head" ${array[$count]} > head2move
			cat "$p" head2move > $p.tail_blunted
			sed -i "1,${lines_head}d" ${array[$count]} # head_blunted
			rm head2move
		fi
		count=$(($count+1))

done < files

blunted_file_count=$(ls *blunted| wc -l | cut -f1 -d ' ')
# echo $blunted_file_count
if [ "$blunted_file_count" = "$file_count" ]; then
	echo "all files blunted correctly"
elif [ "$blunted_file_count" != "$file_count" ]; then
	echo "ERROR: not all files blunted correctly, likely one luck file OK"
fi

# ls x* > allfiles
# while read p; do
# 	echo "$p head tail"
# 	head -3 $p
# 	tail -3 $p
# done < allfiles

mkdir perlrun
mv *blunted perlrun/
cd perlrun
ls x*blunted > bfiles
### run using this index file

parallel -a bfiles perl /work/share/kdy_caolijun/software/popoolation2_1201/fst-sliding.pl --input {} --output {}.ngm90.pairs.sorted.fst --suppress-noninformative --min-count 1 --min-coverage 50 --max-coverage $4 --min-covered-fraction 0 --window-size $6 --step-size $7 --pool-size $5 

### cat the output
#indexname=test
#head -1 xaa*pwc > pwc_header
#head -1 xaa*rc > rc_header
#sed -i '1d' x* 
cat x*params > $outfile.fisher.params
cat x*fst > $outfile.fst
#cat rc_header x*rc > $outfile..allelefreq_rc
rm x*
mv $outfile* ../../
cd ../../
rm -rf temp_popool_fstsliding
