#!/bin/bash

FILENAME=$1
TEMPFILE=mktemp

if [ "$FILENAME" == "" ]; then
	echo "SYNAPSIS: $0 <filename>"
	echo 
	echo "outputs a histogram to <filename>.png"
	exit
fi

{ echo $FILENAME; cat $FILENAME; } > $TEMPFILE

R --slave --quiet --no-save << EOF
png("$FILENAME.png")
data = read.table("$FILENAME", col.names = ("data"))\$data
hist(data, main="$FILENAME")
EOF

echo "histogram is at $FILENAME.png"

rm $TEMPFILE

