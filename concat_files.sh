#!/bin/bash

# usage: sh concat_files.sh <folder> <output file> <sample names>

folder=$1
outfile=$2
names=$3

#5th column(coverage) of the coverage file is extracted from all files and appended to a single file
awk 'BEGIN{FS=OFS="\t"}; { a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i] }' $(ls -1v $folder/*) > $folder/temp.txt
awk '{print $1,"\t",$2,"\t",$3,"\t",$4}' $(ls $folder/* | grep -v $folder/temp.txt |  head -1) > $folder/temp2.txt
paste  $folder/temp2.txt $folder/temp.txt > $folder/test

sort  -V -k1 -k2n -k3n $folder/test > $folder/temp3.txt
echo -e $names > $outfile
cat $folder/temp3.txt >> $outfile


rm $folder/temp.txt $folder/temp2.txt $folder/test $folder/temp3.txt 


