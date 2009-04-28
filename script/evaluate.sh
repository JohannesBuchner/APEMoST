#!/bin/bash
function die() {
	echo "died:" $1
	exit $2
}

TMP=$(mktemp)
mcmcdir=~/Desktop/Arbeit
matrix_man=$mcmcdir/matrix_man.exe


bash $mcmcdir/script/join.sh *-chain*.prob.dump || die 'join' 1
#paste *.all > data || exit 2
for j in *.all; do 
	$mcmcdir/histogram_tool.exe 500 $j | sed 's/\(.*\)\.\..*:[^0-9]*\([0-9]*\)/\1\t\2/g'|grep '\.' > $j.hist.tab || die 'histogram' 3
done
{
echo counting misses 1>&2
for j in Amplitude Phase Frequenz; do
#	cut -f1 $j.all.hist.tab |grep -v $(cut -f1 $mcmcdir/idl/simplesin/correct-$j.all.dat) > $TMP
	cat $j.all|$matrix_man addc i=1 -$(cut -f1 $mcmcdir/idl/simplesin/correct-$j.all.dat) > $TMP

	$mcmcdir/sum_tool.exe -s $TMP   |
		$matrix_man new     i=2 |
		$matrix_man addc    i=2 $(cat $j.all|wc -l)|
		$matrix_man invm    i=2 |
		$matrix_man mul_rel i=1 1|
		$matrix_man rm      i=2
done

echo stuck chains 1>&2
chains=$(cat acceptance_rate.dump |~/Desktop/Arbeit/matrix_man.exe diff 'i>1'|grep -w '0.000000e+00' -c)
echo $(($chains*1000)) # stuck chains are expensive

} | tee $TMP 1>&2

$mcmcdir/sum_tool.exe $TMP

rm $TMP

