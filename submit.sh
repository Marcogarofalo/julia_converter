confs=(`ls /leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE`)
for((ic=0; ic<${#confs[@]}; ic++))
do
	sbatch scr_serial.sh ${confs[ic]}
done
