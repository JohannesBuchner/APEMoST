#!/bin/bash

commands_file=$1
a=$(mktemp)
b=$(mktemp)


cat > $a

while read line; do
	$line < $a > $b
	mv $b $a
done < $commands_file
cat $a
rm $a

