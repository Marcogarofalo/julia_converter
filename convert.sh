contrs_in=("A0P5" "P5P5" "VKVK")
contrs_out=("P5A0" "P5P5" "VKVK")
TMOS=("OS" "TM")
eqop=("equal" "opposite")

ens=("B.72.64" "0.00072" "0.018" "0.020" "mix_fixed"   "end_line"
    "C.06.112" "0.00060" "0.016" "0.018" "mix_fixed"  "end_line"
    "C.06.80" "0.00060" "0.016" "0.018"   ""          "end_line" 
    "D.54.96" "0.00054" "0.013" "0.014"   "mix"          "end_line"  )
#  "B.72.96" "0.00072" "0.018" "0.020"   
#  "E.44.112" "0.00072"   "0.018"  "0.020"
prefixs=("mix_fixed"  )

dir="/leonardo_scratch/large/userexternal/mgarofal/gmtemp"
dirout=$1

tot_ens=`for i in ${ens[@]}; do echo $i ; done | grep end_line  | wc -l`
cols=$((${#ens[@]}/$tot_ens))

for ((ie=0;ie<${#ens[@]};ie+=$cols ))
do

for ((ic=0;ic<${#contrs_in[@]};ic++ ))
do

for ((it = 0; it < ${#eqop[@]}; it++)); do

    #echo a
    #julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/data/${ens[ie+4]}_l_l_${TMOS[it]}_${contrs_in[ic]}    ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_${contrs_out[ic]}.txt
    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/data/${ens[ie+4]}_l_s1_${TMOS[it]}_${contrs_in[ic]}   ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_mu.${ens[ie+2]}_${contrs_out[ic]}.txt
    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/data/${ens[ie+4]}_l_s2_${TMOS[it]}_${contrs_in[ic]}   ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+1]}_mu.${ens[ie+3]}_${contrs_out[ic]}.txt
    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/data/${ens[ie+4]}_s1_s1_${TMOS[it]}_${contrs_in[ic]}  ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+2]}_${contrs_out[ic]}.txt
    julia convert_sanfo_gm2.jl ${dir}/${ens[ie]}/data/${ens[ie+4]}_s2_s2_${TMOS[it]}_${contrs_in[ic]}  ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+3]}_${contrs_out[ic]}.txt

done
done
done
