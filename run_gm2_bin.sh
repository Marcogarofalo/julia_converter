#!/bin/bash

# julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.00072_P5A0.txt  ../flow/data_from_gm2/B96_l1_l1_P5A0.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
# julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.00072_P5P5.txt  ../flow/data_from_gm2/B96_l1_l1_P5P5.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
# julia  convert_gm2_2_bin_no_Conflist.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.0006675_P5A0.txt  ../flow/data_from_gm2/B96_l2_l2_P5A0.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675
# julia  convert_gm2_2_bin_no_Conflist.jl ../g-2_new_stat/cB.72.96_r.opposite_mu.0.0006675_P5P5.txt  ../flow/data_from_gm2/B96_l2_l2_P5P5.dat ../flow/data_from_gm2/B.72.96/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675


julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.64_r.opposite_mu.0.00072_P5A0.txt  ../flow/data_from_gm2/B64_TM_l1_l1_P5A0.dat ../flow/data_from_gm2/B.72.64/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.64_r.opposite_mu.0.00072_P5P5.txt  ../flow/data_from_gm2/B64_TM_l1_l1_P5P5.dat ../flow/data_from_gm2/B.72.64/confsListll.txt  1.77800000000 0.13942650  0.00072 0.00072
julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.64_r.opposite_mu.0.0006675_P5A0.txt  ../flow/data_from_gm2/B64_TM_l2_l2_P5A0.dat ../flow/data_from_gm2/B.72.64/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675
julia  convert_gm2_2_bin.jl ../g-2_new_stat/cB.72.64_r.opposite_mu.0.0006675_P5P5.txt  ../flow/data_from_gm2/B64_TM_l2_l2_P5P5.dat ../flow/data_from_gm2/B.72.64/confsListll.txt  1.77800000000 0.13942650  0.0006675 0.0006675


julia convert_blind_format_2bin.jl ../gm2-mistuning/P5P5_LMA_confs/B64/P5P5_tm_mu7.2000e-04_LMA.bin \
        ../flow/data_fpi_A0/B64_TM_l1_l1_P5P5.dat  ../gm2-mistuning/P5P5_LMA_confs/B64/B64_ncfg1382_neta1024.txt  \
        1.77800000000 0.13942650  0.00072 0.00072
julia convert_blind_format_2bin.jl ../gm2-mistuning/valence_mu_light/B64/P5P5_tm_mistuning_B64.bin \
        ../flow/data_fpi_A0/B64_TM_deriv_P5P5.dat  ../gm2-mistuning/valence_mu_light/B64/B64_mu_val_conf.txt \
        1.77800000000 0.13942650  0.0006675 0.0006675



julia convert_blind_format_2bin.jl ../gm2-mistuning/P5P5_LMA_confs/B64/P5P5_OS_mu7.2000e-04_LMA.bin \
        ../flow/data_fpi_A0/B64_OS_l1_l1_P5P5.dat  ../gm2-mistuning/P5P5_LMA_confs/B64/B64_ncfg1382_neta1024.txt  \
        1.77800000000 0.13942650  0.00072 0.00072
julia convert_blind_format_2bin.jl ../gm2-mistuning/valence_mu_light/B64/P5P5_OS_mistuning_B64.bin \
        ../flow/data_fpi_A0/B64_OS_deriv_P5P5.dat  ../gm2-mistuning/valence_mu_light/B64/B64_mu_val_conf.txt \
        1.77800000000 0.13942650  0.0006675 0.0006675


#### B24
julia convert_small_vol_2bin.jl ../flow/data_fpi_A0/smallV/B24/A0P5_tm_mu2.5000e-03_mu2.5000e-03.bin  ../flow/data_fpi_A0/B24_TM_l1_l1_P5A0.dat /dev/null 1.778 0.1394267 0 0
julia convert_small_vol_2bin.jl ../gm2-mistuning/B24/P5P5_tm_mu2.5000e-03_mu2.5000e-03.bin            ../flow/data_fpi_A0/B24_TM_l1_l1_P5P5.dat /dev/null 1.778 0.1394267 0 0