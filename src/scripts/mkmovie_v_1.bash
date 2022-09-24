#!/bin/bash
#@ mkmovie.bash::
#@ Eraldo Pereira Marinho
#@ Mon Sep 17 18:54:53 BRT 2012
#@ Version 1
#@ Revision 1.2
#@ Revision 1.3
objname="$1"
#data frame is defined as the rectangle
# x in [a:b], and y in [c:d]
a="$2"
b="$3"
c="$4"
d="$5"
# x-coordinate column position in the data file
x="$6"
# y-coordinate column position in the data file
y="$7"
# number of horizontal pixels in the output image
sx="$8"
# resize geometry
Sx="${9}"
# initial frame number
m="${10}"
# final frame number
n="${11}"
# color modifier (leave empty if you don't know what does it mean)
if [ -z ${12} ]
then
	fps=25
	dcolor=""
else
	fps=${12}
	dcolor="${13}"
fi
#
plotxy=~/workspace/programs/rasplot/bin/plotxy$dcolor
#
let n++
for((i=m;i<n;i++))
do
	let j=100000+i
	datafile=$objname$i.data
	pngfile=$objname$j.png
	if [ ! -e $pngfile ]
	then
		if [ -e $datafile.gz ]
		then
			gunzip $datafile.gz
			gzflag=1
		else
			gzflag=
		fi
		if [ -e $datafile -a -s $datafile ]
		then
			echo $datafile '->' $pngfile
			rasterfile=$objname$j.ras
			$plotxy $datafile $a $b $c $d $x $y $sx > $rasterfile
			[ -z $gzflag ] || gzip -f $datafile
			convert -resize ${Sx}x${Sx} $rasterfile $pngfile
			rm -f $rasterfile
		fi
	fi
done
mencoder "mf://$objname*.png" -noskip -mf fps=$fps -o $objname.avi -ovc lavc -lavcopts vcodec=mpeg4:autoaspect=1
