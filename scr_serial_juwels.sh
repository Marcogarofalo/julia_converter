#!/bin/sh
#SBATCH -A libegm2
# #SBATCH -p lrd_all_serial
#SBATCH --time 04:00:00     # format: HH:MM:SS

##SBATCH -p boost_usr_prod
##SBATCH --qos boost_qos_dbg
##SBATCH --time 00:30:00     # format: HH:MM:SS

#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node=1 # 4 tasks out of 32
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=1 # 4 tasks out of 32
#SBATCH --job-name=convertLIBE
#SBATCH --mem=5000MB          # memory per node out of 494000MB (481GB)

# julia   convert_one_LIBE.jl  $1 $2
julia convert_one_LIBE.jl 0008_r1   input_LIBE_B64.jl
# confs=(`ls /leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE`)
# for((ic==0; ic<${#confs[@]}; ic+=4))
# do
#     julia   convert_one_LIBE.jl ${confs[ic]}   &
#     pid0=$!
#     julia   convert_one_LIBE.jl ${confs[ic+1]} &
#     pid1=$!
#     julia   convert_one_LIBE.jl ${confs[ic+2]} &
#     pid2=$!
#     julia   convert_one_LIBE.jl ${confs[ic+3]}  &
#     pid3=$!
#     wait $pid0 $pid1 $pid2 $pid3
# done 

