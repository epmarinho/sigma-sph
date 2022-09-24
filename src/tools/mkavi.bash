#!/bin/bash
objname="$1"
m=$2
n=$3
let n++
for((i=m;i<n;i++))
do
	let j=10000+i
	mv $objname$i.png out$j.png
done
mencoder "mf://out*.png" -mf fps=30 -o $objname.avi -ovc lavc -lavcopts vcodec=mpeg4
