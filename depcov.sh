#!/bin/bash

cd $1
touch total.txt

for dir in $(ls)
do
if [ -d "$dir" ]
then
	ID=$dir
	cd $dir
	cd mapping_data
	depth_tmp=$(samtools depth $ID"_mapped_sorted_rmdup.bam" |  awk '{sum+=$3} END { print sum}')
	depth=$(bc<<<"scale=3;$depth_tmp/3099734149")
	
	coverage=$(bedtools bamtobed -i $ID"_mapped_sorted_rmdup.bam" | awk '{sum+=$3-$2} END {print sum/3099734149*100}')
	echo $ID "Depth" $depth"X" "Coverage" $coverage"%" >> ../../total.txt
	cd ../../
fi
done
