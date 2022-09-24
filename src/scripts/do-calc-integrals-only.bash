#!/bin/bash
# Mon Jun 20 11:31:12 AM -03 2022
# Use this script to compute integral quantities outcoming from
# a series of SPH simulations
simcode=./bin/sigma-treecode_calc_integrals_only
seriesname=$1
# initial time-level
m=$2
# final time-level
n=$3
# number K of nearest neighbors
K=$4
# time-step
dt=$5
# squared aperture parameter
theta2=$6
# squared softening-length
epsilon2=$7

if [ -e $seriesname$m.data.gz ]
then
	gunzip $seriesname$m.data.gz
	gzflag=1
fi
for((i=m;i<n;i++))
do
	input_file=$seriesname$i.data
	[ -s $input_file ] && "$simcode" $input_file $K $dt $theta2 $epsilon2
	[ -z $gzflag ] || gzip $seriesname$i.data
done
[ -z $gzflag ] || gzip $seriesname$m.data
