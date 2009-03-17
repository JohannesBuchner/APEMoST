#!/bin/bash

FILENAME=$1
TEMPFILE=mktemp

{ echo $FILENAME; cat $FILENAME; } > $TEMPFILE

R --slave --quiet --no-save << EOF
png("$FILENAME.png")
data = read.table("$FILENAME", col.names = ("data"))\$data
hist(data, main="$FILENAME")
EOF

echo "histogram is at $FILENAME.png"

rm $TEMPFILE

