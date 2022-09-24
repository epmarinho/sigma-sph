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
# number of horizontal pixels in the output image used in plotxy
sx="$8"
# convert resize geometry parameter
Sx="${9}"
# color modifier (leave empty if you don't know what does it mean)
if [ -z ${10} ]
then
	fps=25
	dcolor=""
else
	dcolor="${10}"
	colorcolumn="${11}"
	fps=${12}
fi
plotxy=~/workspace/programs/rasplot/bin/plotxy$dcolor
echo $plotxy
for datafile in $objname*.data
do
	pngfile=`basename $datafile data`png
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
			rasterfile=`basename $datafile data`ras
			$plotxy $datafile $a $b $c $d $x $y $colorcolumn $sx > $rasterfile
			[ -z $gzflag ] || gzip -f $datafile
			convert -resize ${Sx}x${Sx} $rasterfile $pngfile
			rm -f $rasterfile
		fi
	fi
done
mencoder "mf://$objname%d.png" -noskip -mf fps=$fps -o $objname.avi -ovc lavc -lavcopts vcodec=mpeg4:autoaspect=1
