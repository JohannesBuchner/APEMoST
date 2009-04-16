
paramnames="Amplitude Phase Frequenz"

function getlimits() {
	case "$1" in
		Amplitude.all)
			echo "[0:3]"
			;;
		Frequenz.all)
			#echo "[14.8:15.9]"
			echo "[14:16]"
			;;
		Phase.all)
			echo "[0:7]"
			;;
	esac
}

mcmcdir=~/Desktop/Arbeit/

for dir in *_*; do 
	cd $dir || continue
	echo $dir
	bash $mcmcdir/script/join.sh *.dump || exit 1
	#paste *.all > data || exit 2
	for j in *.all; do 
		$mcmcdir/histogram_tool.exe 500 $j | sed 's/\(.*\)\.\..*:[^0-9]*\([0-9]*\)/\1\t\2/g'|grep '\.' > $j.hist.tab || exit 3
		rm $j
		echo 'set title "'$dir-$j'";set terminal png; set output "'$j.hist.tab.png'"; set term png size 400, 200; ' plot "$(getlimits $j) "'"'$j.hist.tab'" notitle with boxes, "'$mcmcdir/idl/simplesin/correct-$j.dat'" notitle linecolor 3 with points;'|gnuplot || exit 4
		rm $j.hist.tab
		pngtopnm $j.hist.tab.png > ../${dir}_${j}.hist.tab.pnm || exit 5
		rm $j.hist.tab.png
	done
	cd - > /dev/null
done

echo joining graphs

for name in $paramnames
do
	#*_*_${name}.all.hist.tab.pnm 
	pnmcat -topbottom {{equidistant,chebyshev}_beta,{equidistant,chebyshev}_temperature}_${name}.all.hist.tab.pnm > $name.topbottom.pnm || exit 6
	rm *_*${name}.all.hist.tab.pnm
done

pnmcat -leftright *.topbottom.pnm | pnmtopng > leftright.png || exit 7
rm *.topbottom.pnm
echo "$PWD/leftright.png done"
