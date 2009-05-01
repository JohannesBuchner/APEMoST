#!/bin/bash

cat="cat"
if [ "$1" == "--last" ]; then
	shift
	cat="tail -n $1"
	shift
fi
if [ "$1" == "--first" ]; then
	shift
	cat="head -n $1"
	shift
fi

FILENAME=$1
TEMPFILE=mktemp

if [ "$FILENAME" == "" ]; then
	echo "SYNAPSIS: $0 [--last <n>|--first <n>] <filename>"
	echo 
	echo "outputs a histogram to <filename>.png"
	exit
fi

$cat $FILENAME > $TEMPFILE

echo $($cat $FILENAME|wc -l) datapoints

R --slave --quiet --no-save << EOF
png("$FILENAME.png")
data = read.table("$TEMPFILE", col.names = ("data"))\$data
hist(data, main="$FILENAME")
EOF

echo "histogram is at $FILENAME.png"

rm $TEMPFILE

