#!/bin/bash
function die() {
	echo "died:" $1
	exit $2
}

mcmcdir=~/Desktop/Arbeit

for dir in $*
do
	echo $dir
	TMP=$(mktemp)
	for par in Amplitude Phase Frequenz
	do 
		file=$dir/$par*chain0*.dump
		min=$(grep $par $mcmcdir/idl/simplesin/params|cut -f2)
		max=$(grep $par $mcmcdir/idl/simplesin/params|cut -f3)
		real=$(cut -f1  $mcmcdir/idl/simplesin/correct-$par*)
		#echo -n $par"       "
		$mcmcdir/peaks.exe $min $max $file|head -n2|tail -n1|
			awk '{print (($1-'$real')/('$max' - '$min'))^2*100000"\t"($2+$3)/('$max' - '$min') * 1000"\t"(1-$4)^4*10000;}'|tee -a $TMP
	done
	$mcmcdir/sum_tool.exe -1 $TMP
	rm $TMP

done

