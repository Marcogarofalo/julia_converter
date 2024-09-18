contrs_in=("A0P5" "P5P5" "VKVK")
contrs_out=("P5A0" "P5P5" "VKVK")
TMOS=("OS" "TM")
eqop=("equal" "opposite")

ens=(
    # "C.06.112"  "0.00060"   "0.013"  "0.014"  "0.18000" "0.19000" "0.20000"      "1.6067e-02"  "C112lsc/ave/data"  "end_line"
    # "D.54.96"  "0.00054"   "0.013"  "0.014"  "0.15000" "0.16000" "0.17000"       "1.3557e-02"  "D96lsc/ave/data"  "end_line"
    # "E.44.112" "0.00044"   "0.011"  "0.012"  "0.13000" "0.14000" "0.15000"       "1.1759e-02"  "E112lsc/ave/data"  "end_line"
)
prefixs=("mix_fixed")

# dir="/leonardo_scratch/large/userexternal/mgarofal/gmtemp"
dir="/leonardo_scratch/large/userexternal/fsanfili/Vus/"
dirout=$1

tot_ens=$(for i in ${ens[@]}; do echo $i; done | grep end_line | wc -l)
cols=$((${#ens[@]} / $tot_ens))

for ((ie = 0; ie < ${#ens[@]}; ie += $cols)); do

    for ((ic = 0; ic < ${#contrs_in[@]}; ic++)); do

        for ((it = 0; it < ${#eqop[@]}; it++)); do

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_c1_c1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 4]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_c2_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_c3_c3_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 6]}_${contrs_out[ic]}.txt

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_s_c1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 7]}_mu.${ens[ie + 4]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_s_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 7]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 8]}/_s_c3_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 7]}_mu.${ens[ie + 6]}_${contrs_out[ic]}.txt

        done
    done
done

########################### C80
ens=(
    # "C.06.80" "0.00060" "0.0005850" "0.014" "0.18000" "0.19000" "there is no 3rd charm" "1.6000e-02" "1.7000e-02" "C80lsc/ave/data" "end_line"
      "B.72.64" "0.00720" "0.0006675" "0.014" "0.23000" "0.24000" "there is no 3rd charm" "1.7500e-02" "1.8500e-02" "B64lsc/ave/data" "end_line"


)
prefixs=("mix_fixed")

# dir="/leonardo_scratch/large/userexternal/mgarofal/gmtemp"
dir="/leonardo_scratch/large/userexternal/fsanfili/Vus/"
dirout=$1

tot_ens=$(for i in ${ens[@]}; do echo $i; done | grep end_line | wc -l)
cols=$((${#ens[@]} / $tot_ens))

for ((ie = 0; ie < ${#ens[@]}; ie += $cols)); do

    for ((ic = 0; ic < ${#contrs_in[@]}; ic++)); do

        for ((it = 0; it < ${#eqop[@]}; it++)); do

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_c1_c1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 4]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_c2_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt
            # julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_c3_c3_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 6]}_${contrs_out[ic]}.txt

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_s1_c1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 7]}_mu.${ens[ie + 4]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_s1_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 7]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_s2_c1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 8]}_mu.${ens[ie + 4]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_s2_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 8]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt

            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_l1_s1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 1]}_mu.${ens[ie + 7]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_l1_s2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 1]}_mu.${ens[ie + 8]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_l2_s1_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 2]}_mu.${ens[ie + 7]}_${contrs_out[ic]}.txt
            julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_l2_s2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 2]}_mu.${ens[ie + 8]}_${contrs_out[ic]}.txt

            # julia convert_sanfo_gm2.jl ${dir}/${ens[ie + 9]}/_s2_c2_${TMOS[$(((it)))]}_${contrs_in[ic]} ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie + 8]}_mu.${ens[ie + 5]}_${contrs_out[ic]}.txt

            # julia convert_sanfo_gm2.jl ${dir}/${ens[ie+8]}/_s_c3_${TMOS[$(((it)))]}_${contrs_in[ic]}  ${dirout}/c${ens[ie]}_r.${eqop[it]}_mu.${ens[ie+7]}_mu.${ens[ie+6]}_${contrs_out[ic]}.txt

        done
    done
done

# # create fake charm
# for ((ic=0;ic<${#contrs_in[@]};ic++ )); do
# for ((it = 0; it < ${#eqop[@]}; it++)); do
#     julia create_fake_gm2.jl 224  ${dirout}/cC.06.112_r.${eqop[it]}_mu.0.99999_${contrs_out[ic]}.txt
#     julia create_fake_gm2.jl 224  ${dirout}/cE.44.112_r.${eqop[it]}_mu.0.99999_${contrs_out[ic]}.txt
# done
# done
