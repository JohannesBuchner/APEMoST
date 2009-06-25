#!/bin/bash
function die() {
	echo "died:" $1
	exit $2
}

mcmcdir=~/Desktop/Arbeit
measure=$mcmcdir/peaks.exe
if [ "$CATCOMMAND" == "" ]; then
	cat="cat"
else
	cat="$CATCOMMAND"
fi

TMP=$(mktemp)
cat params|cut -f4|while read param; do 
	echo -n $param "     "; 
	min=$(grep -w "$param" params|cut -f2); 
	max=$(grep -w "$param" params|cut -f3); 
	$cat $param-chain-0.prob.dump > $TMP
	$mcmcdir/peaks.exe $min $max $TMP|head -n2|tail -n1
done

