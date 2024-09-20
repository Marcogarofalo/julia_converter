contrs_in=("A0P5" "P5P5" "VKVK")
contrs_out=("P5A0" "P5P5" "VKVK")
TMOS=("OS" "TM")
eqop=("equal" "opposite")

ens=(
    # "B.72.64" "0.0006675" "0.00072" "0.018254"  "mix"   "end_line"
    # "C.06.80" "0.000585" "0.00060"  "0.016067"  "mix"   "end_line" 
    # "D.54.96" "0.0004964" "0.00054" "0.013557"  "mix"   "end_line"  
    "E.44.112" "4.3100e-4" "0.00044"  "1.1759e-02"   "mix"   "end_line"     
     ### "C.06.112" "0.00060" "0.016" "0.018" "mix_fixed"  "end_line"
    # "B.72.96" "0.00072" "0.018" "0.019"  ""           "end_line"   
     )


# dir="/leonardo_scratch/large/userexternal/mgarofal/gmtemp"
dir="/leonardo_work/INF24_lqcd123_1/mgarofal/tau_decay_strange_bin_mu_corr"
dirout=$1

tot_ens=`for i in ${ens[@]}; do echo $i ; done | grep end_line  | wc -l`
cols=$((${#ens[@]}/$tot_ens))

for ((ie=0;ie<${#ens[@]};ie+=$cols ))
do

for ((ic=0;ic<${#contrs_in[@]};ic++ ))
do

for ((it = 0; it < ${#eqop[@]}; it++)); do

    #echo a
    # if [ ${ens[ie]} == "E.44.112" ]; then
    #     julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/${ens[ie+4]}_l_l_${TMOS[it]}_${contrs_in[ic]}    ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_${contrs_out[ic]}.txt
    # else
    #     julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/ll_${TMOS[it]}_${contrs_in[ic]}    ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_${contrs_out[ic]}.txt
    # fi

    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/${ens[ie+4]}_l1_s_${TMOS[it]}_${contrs_in[ic]}   ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_mu.${ens[ie+3]}_${contrs_out[ic]}.txt
    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/${ens[ie+4]}_l2_s_${TMOS[it]}_${contrs_in[ic]}   ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+2]}_mu.${ens[ie+3]}_${contrs_out[ic]}.txt
    # julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/${ens[ie+4]}_s1_s1_${TMOS[it]}_${contrs_in[ic]}  ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+2]}_${contrs_out[ic]}.txt
    # julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/${ens[ie+4]}_s2_s2_${TMOS[it]}_${contrs_in[ic]}  ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+3]}_${contrs_out[ic]}.txt

    
done
done
done


# create fake charm
# for ((ic=0;ic<${#contrs_in[@]};ic++ )); do
# for ((it = 0; it < ${#eqop[@]}; it++)); do
#     julia create_fake_gm2.jl 224  ${dirout}/cC.06.112_r.${eqop[it]}_mu.0.99999_${contrs_out[ic]}.txt
#     julia create_fake_gm2.jl 224  ${dirout}/cE.44.112_r.${eqop[it]}_mu.0.99999_${contrs_out[ic]}.txt
# done
# done