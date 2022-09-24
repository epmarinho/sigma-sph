#!/bin/bash
simcode=$1
seriesname=$2
# initial time-level
m=$3
# final time-level
n=$4
# number K of nearest neighbors
K=$5
# time-step
dt=$6
# squared aperture parameter
theta=$7
# squared softening-length
epsilon=$8
# number of OMP threads
numthreads=$9

if [ -e $seriesname$m.data.gz ]
then
	gunzip $seriesname$m.data.gz
	gzflag=1
fi
for((i=m;i<n;i++))
do
	let j=i+1
	#echo $i '->' $j
	input_file=$seriesname$i.data
	[ -s $input_file ] && time "$simcode" $input_file $K $dt $theta $epsilon $numthreads  >$seriesname$j.data
	[ -z $gzflag ] || gzip $seriesname$i.data
done
[ -z $gzflag ] || gzip $seriesname$m.data

# Last change:
# dom 12 dez 2021 10:16:07 -03 @capricornio.rc.unesp.br
