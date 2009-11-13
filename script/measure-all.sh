#!/bin/bash

# 
# @author Johannes Buchner
# This script, run in a working directory with a params file and
# dumpfiles, evaluates the dump files using peaks.exe for every 
# parameter in params.
# 
# called without parameters.
# 



mcmcdir=$(dirname $0)/..
measure=$mcmcdir/peaks.exe
if [ "$CATCOMMAND" == "" ]; then
	cat="cat"
else
	cat="$CATCOMMAND"
fi

TMP=$(mktemp)
cat params|awk '{ print $4; }'|while read param; do 
	echo $param "     "; 
	min=$(grep -w "$param" params|awk '{ print $2; }'); 
	max=$(grep -w "$param" params|awk '{ print $3; }'); 
	$cat $param-chain-0.prob.dump | head -n-1 > $TMP
	$mcmcdir/peaks.exe $min $max $TMP|{ read; cat; }
done
rm $TMP

