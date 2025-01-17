#!/bin/bash

# confs=(`ls /leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE`)
confs=(`ls /p/scratch/libegm2/bacchio1/B64/pion_LIBE/`)
for((ic=0; ic<${#confs[@]}; ic++))
do
	sbatch scr_serial_juwels.sh ${confs[ic]} input_LIBE_B64.jl
done
