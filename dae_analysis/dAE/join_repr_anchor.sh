#!/bin/bash

FN=$1

awk -F "\t" 'BEGIN {
while (getline < "'"$FN"'"){
	A[$1] = $2;
}
close("'"$FN"'")
}
{
	split($4, cell_tf, "_");
	if($4 in A)
		print $0"\t"A[$4]
}' $2

