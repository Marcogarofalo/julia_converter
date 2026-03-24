#!/bin/bash

# julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.00072_P5A0.txt  ../flow/data_from_gm2/B96_l1_l1_P5A0.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
# julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.00072_P5P5.txt  ../flow/data_from_gm2/B96_l1_l1_P5P5.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
# julia  convert_gm2_2_bin_no_Conflist.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.0006675_P5A0.txt  ../flow/data_from_gm2/B96_l2_l2_P5A0.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675
# julia  convert_gm2_2_bin_no_Conflist.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.0006675_P5P5.txt  ../flow/data_from_gm2/B96_l2_l2_P5P5.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675

mus=(0.00072  0.0006675)
reg=("tm"  "OS")
reg_in=("opposite"  "equal")
gammas=("A0P5"  "P5P5")
gammas_in=("P5A0"  "P5P5")

for ((imu=0; imu<${#mus[@]}; imu++)); do
for ((ir=0; ir<${#reg[@]}; ir++)); do
for ((ig=0; ig<${#gammas[@]}; ig++)); do
julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.64_r.${reg_in[$ir]}_mu.${mus[$imu]}_${gammas_in[$ig]}.txt  ../flow/data_fpi_A0/B64_${reg[$ir]}_l${imu}_l${imu}_${gammas[$ig]}.dat ../flow/data_from_gm2/B.72.64/confsListll.txt  1.77800000000 0.13942650  ${mus[$imu]} ${mus[$imu]}
done done done

mus=(0.00060 0.0005850)
for ((imu=0; imu<${#mus[@]}; imu++)); do
for ((ir=0; ir<${#reg[@]}; ir++)); do
for ((ig=0; ig<${#gammas[@]}; ig++)); do
# julia  convert_gm2_2_bin.jl ../flow/data_from_gm2/C.06.80/cC.06.80_r.${reg_in[$ir]}_mu.${mus[$imu]}_${gammas_in[$ig]}.txt  ../flow/data_fpi_A0/C80_${reg[$ir]}_l${imu}_l${imu}_${gammas[$ig]}.dat ../flow/data_from_gm2/C.06.80/confsListMix  1.83600000000 0.138752850000  ${mus[$imu]} ${mus[$imu]}
julia  convert_gm2_2_bin.jl ../g-2_new_stat//cC.06.80_r.${reg_in[$ir]}_mu.${mus[$imu]}_${gammas_in[$ig]}.txt  ../flow/data_fpi_A0/C80_${reg[$ir]}_l${imu}_l${imu}_${gammas[$ig]}.dat ../flow/data_from_gm2/C.06.80/confsList  1.83600000000 0.138752850000  ${mus[$imu]} ${mus[$imu]}

done done done

# mus=(0.00060  0.0005850)
# for ((imu=0; imu<${#mus[@]}; imu++)); do
# for ((ir=0; ir<${#reg[@]}; ir++)); do
# for ((ig=0; ig<${#gammas[@]}; ig++)); do
# julia  convert_gm2_2_bin.jl ../g-2_new_stat/cC.60.80_r.${reg_in[$ir]}_mu.${mus[$imu]}_${gammas_in[$ig]}.txt  ../flow/data_fpi_A0/C80_${reg[$ir]}_l${imu}_l${imu}_${gammas[$ig]}.dat ../flow/data_from_gm2/C.60.80/confsListll.txt  1.77800000000 0.13942650  ${mus[$imu]} ${mus[$imu]}
# done done done


# julia convert_blind_format_2bin.jl ../gm2-mistuning/P5P5_LMA_confs/B64/P5P5_tm_mu7.2000e-04_LMA.bin \
#         ../flow/data_fpi_A0/B64_TM_l0_l0_P5P5.dat  ../gm2-mistuning/P5P5_LMA_confs/B64/B64_ncfg1382_neta1024.txt  \
#         1.77800000000 0.13942650  0.00072 0.00072
# julia convert_blind_format_2bin.jl ../gm2-mistuning/valence_mu_light/B64/P5P5_tm_mistuning_B64.bin \
#          ../flow/data_fpi_A0/B64_tm_val_deriv_P5P5.dat  ../gm2-mistuning/valence_mu_light/B64/B64_mu_val_conf.txt \
#          1.77800000000 0.13942650  0.00072 0.0006675



# julia convert_blind_format_2bin.jl ../gm2-mistuning/P5P5_LMA_confs/B64/P5P5_OS_mu7.2000e-04_LMA.bin \
#         ../flow/data_fpi_A0/B64_OS_l1_l1_P5P5.dat  ../gm2-mistuning/P5P5_LMA_confs/B64/B64_ncfg1382_neta1024.txt  \
#         1.77800000000 0.13942650  0.00072 0.00072
# julia convert_blind_format_2bin.jl ../gm2-mistuning/valence_mu_light/B64/P5P5_OS_mistuning_B64.bin \
#         ../flow/data_fpi_A0/B64_OS_deriv_P5P5.dat  ../gm2-mistuning/valence_mu_light/B64/B64_mu_val_conf.txt \
#         1.77800000000 0.13942650  0.00072 0.0006675


#### B24
gammas="A0P5  P5P5"
# for g in $gammas; do
mus=(2.5000e-03  4.2000e-03  5.9000e-03)
reg="tm  OS"
for g in $gammas; do
for ((i=0; i<${#mus[@]}; i++)); do
for r in $reg; do
julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/B24/${g}_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/B24_${r}_l${i}_l${i}_${g}.dat /dev/null 1.778 0.1394267 0 0
#julia convert_small_vol_2bin.jl ../gm2-mistuning/B24/P5P5_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/B24_${r}_l${i}_l${i}_P5P5.dat /dev/null 1.778 0.1394267 0 0
done done done

# B32
mus=(2.5000e-03  4.2000e-03  5.9000e-03)
for g in $gammas; do
for ((i=0; i<${#mus[@]}; i++)); do
for r in $reg; do
julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/B32/${g}_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/B32_${r}_l${i}_l${i}_${g}.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../gm2-mistuning/B32/${g}_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/B32_${r}_l${i}_l${i}_${g}.dat /dev/null 1.778 0.1394267 0 0
done done done

#C48
mus=(2.0000e-03  3.4000e-03  4.8000e-03)
for g in $gammas; do
for ((i=0; i<${#mus[@]}; i++)); do
for r in $reg; do
julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/C48/${g}_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/C48_${r}_l${i}_l${i}_${g}.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../gm2-mistuning/B32/${g}_${r}_mu${mus[$i]}_mu${mus[$i]}.bin  ../flow/data_fpi_A0/B32_${r}_l${i}_l${i}_${g}.dat /dev/null 1.778 0.1394267 0 0
done done done
# done
# julia convert_small_vol_2bin.jl ../gm2-mistuning/B24/P5P5_tm_mu2.5000e-03_mu2.5000e-03.bin            ../flow/data_fpi_A0/B24_TM_l1_l1_P5P5.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/B24/A0P5_tm_mu4.2000e-03_mu4.2000e-03.bin  ../flow/data_fpi_A0/B24_TM_l2_l2_P5A0.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../gm2-mistuning/B24/P5P5_tm_mu4.2000e-03_mu4.2000e-03.bin            ../flow/data_fpi_A0/B24_TM_l2_l2_P5P5.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/B24/A0P5_tm_mu5.9000e-03_mu5.9000e-03.bin  ../flow/data_fpi_A0/B24_TM_l3_l3_P5A0.dat /dev/null 1.778 0.1394267 0 0
# julia convert_small_vol_2bin.jl ../gm2-mistuning/B24/P5P5_tm_mu5.9000e-03_mu5.9000e-03.bin            ../flow/data_fpi_A0/B24_TM_l3_l3_P5P5.dat /dev/null 1.778 0.1394267 0 0
