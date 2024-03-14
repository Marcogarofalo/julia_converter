using HDF5
using Printf
import Base.write
using Profile
using Statistics

include("./read_hdf5.jl")
# qui i correlatori hanno 5 componenti
# le combinazioni sono 00, 01, 10, 02, 20
# dove 0 e' il propagator normale, 1 e' D^2 e 2 e' D g5 D
include("./gamma.jl")


function hits_SIB(outfile::IOStream,
	confname::String, T::Int32, hits_qcd::Vector{String},
	masses::Array{Float64, 1}, TMOSs::Vector{Vector{String}},
	info_counterterms::Array{Int32, 1}, counterterms::Array{Float64, 1},
	gammas::Vector{String})

	# h5ls -r file gives: {96, 1, 3, 3, 16, 2}
	raw_data::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, 5, 1, T)
	m1list::Array{Float64, 1} = [masses[1], masses[1]]
	corr::Array{Float64, 8} = zeros(Float64, length(hits_qcd), length(masses), 2, length(TMOSs), 16, length(counterterms), T, 2)

	for (ihit, hit) in enumerate(hits_qcd)
		fname::String = string(confname, hit)
		fid = h5open(fname, "r")
		for (im, m) in enumerate(masses)
			m1list[1] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)

					combo::String = @sprintf("%s/mesons/%c%.4e_%c%.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					# raw_data = fid[combo][:, id, 1, :]
					raw_data = read(fid, combo)
					# raw_data = h5read(fname, combo)

					for ig in 1:16

						for ic in 1:length(counterterms)

							for t in 1:T
								corr[ihit, im, im1, iTMOS, ig, ic, t, 1] += raw_data[1, ig, ic, 1, t]
								corr[ihit, im, im1, iTMOS, ig, ic, t, 2] += raw_data[2, ig, ic, 1, t]
							end
							# corr[im, im1, iTMOS, ig, ie, :, 1] .+= raw_data[1, id, e1, e2, 1, :]

						end
					end

				end
			end
		end
		close(fid)
	end

	# corr_qcd ./= length(hits_qcd)

	# writing
	# for reim in 1:2
	# 	for t in 1:T
	# 		for i in 1:length(counterterms)
	# 			for (ig, g) in enumerate(gammas)
	# 				for (iTMOS, TMOS) in enumerate(TMOSs)
	# 					for (im,m) in enumerate(masses)
	# 						for im1 in 1:2
	# 							corr[ihit, im, im1, iTMOS, ig, i, t, reim] /= length(hits_qcd)
	# 							# corr[im, im1, iTMOS, ig, i, t, 2] /= length(hits_qcd)
	# 						end
	# 					end
	# 				end
	# 			end
	# 		end
	# 	end
	# end
	return (corr)
end

function ave_hits_loop(outfile::IOStream,
	confname::String, T::Int32, hits_qcd::Vector{String},
	masses::Array{Float64, 1}, TMOSs::Vector{Vector{String}},
	info_counterterms::Array{Int32, 1}, counterterms::Array{Float64, 1},
	gammas::Vector{String})

	# h5ls -r file gives: {1, 1, 3, 16, 2}
	raw_data::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, 3, 1, 1)
	m1list::Array{Float64, 1} = [masses[1], masses[1]]
	corr::Array{Float64, 7} = zeros(Float64, length(hits_qcd), length(masses), 2, length(TMOSs), 16, length(counterterms), 2)

	for (ihit, hit) in enumerate(hits_qcd)
		fname::String = string(confname, hit)
		fid = h5open(fname, "r")
		for (im, m) in enumerate(masses)
			m1list[1] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)

					combo::String = @sprintf("%s/mesons/%c%.4e_%c%.4e_loop", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					# raw_data = fid[combo][:, id, 1, :]
					raw_data = read(fid, combo)
					# raw_data = h5read(fname, combo)

					for ig in 1:16
						for ic in 1:length(counterterms)


							corr[ihit, im, im1, iTMOS, ig, ic, 1] += raw_data[1, ig, ic, 1, 1]
							corr[ihit, im, im1, iTMOS, ig, ic, 2] += raw_data[2, ig, ic, 1, 1]


						end
					end

				end
			end
		end
		close(fid)
	end


	return (corr)
end




function Gamma_contraction(acc_data::Array{Float64, 5}, raw_data::Array{Float64, 8}, T::Int32, iTMOS::Int)
	for t in 1:T
		
		for ins in 1:5
			for i in 1:16
				for j in 1:4
					for k in 1:4
						tmp::Complex{Float64} = (raw_data[1, j, Gamma[i].col[k], k, Gamma[i].col[j], ins, 1, t]
												 +
												 (raw_data[2, j, Gamma[i].col[k], k, Gamma[i].col[j], ins, 1, t])im) *
												Gamma[i].val[j] * Gamma[i].val[k]
						acc_data[t, i, ins, iTMOS, 1] += tmp.re
						acc_data[t, i, ins, iTMOS, 2] += tmp.im
					end
				end
			end
		end

	end
end

function read_open(conf::String, acc_data::Array{Float64, 5}, basename::String, L::Int32, T::Int32, m::Float64, TMOSs::Vector{Vector{String}})
	# {128, 1, 4, 4, 4, 4, 2}
	raw_data::Array{Float64, 8} = Array{Float64, 8}(undef, 2, 4, 4, 4, 4, 5, 1, T)


	hits = readdir(string(basename, "/", conf))
	pattern = string("^twop2_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new[1:100]]
	println("Nhits:   ", length(hits))
	# println(hits)
	for (i, hit) in enumerate(hits)
		filename::String = string(basename, "/", conf, "/", hit)
		fid::HDF5.File = h5open(filename, "r")
		# println(keys(fid)[1],typeof(fid))
		# println(hit)
		for (iTMOS, TMOS) in enumerate(TMOSs)
			group::String = @sprintf("/%s/mesons/%c%.4e_%c%.4e_open", keys(fid)[1], TMOS[1], m, TMOS[2], m)

			raw_data = read(fid, group)
			Gamma_contraction(acc_data, raw_data, T, iTMOS)


		end
		close(fid)
	end

	factor::Float64 = L^3 * length(hits)
	for (iTMOS, TMOS) in enumerate(TMOSs)
		for t in 1:T
			for i in 1:16
				for ins in 1:5

					acc_data[t, i, ins, iTMOS, 1] /= factor
					acc_data[t, i, ins, iTMOS, 2] /= factor
				end
			end
		end
	end
	# symm_t!(acc_data)

	# println("stoc stoc")
	# for i in 1:16
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 0, acc_data[1, i, 1, 1], acc_data[1, i, 1, 2], acc_data[1, i, 2, 1], acc_data[1, i, 2, 2])
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 1, acc_data[2, i, 1, 1], acc_data[2, i, 1, 2], acc_data[2, i, 2, 1], acc_data[2, i, 2, 2])
	# end
end

function main()
	# basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_SIB"
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/cold"

	T::Int32 = 8
	L::Int32 = 4
	V::Int = L * L * L
	beta = 1.778000000000 ##check
	kappa = 0.14 ##check
	masses = [0.001]

	info_counterterms::Vector{Int32} = [5]
	# e::Float64 = 1.0000e-03 # before there was written by mistake 1e-2  
	counterterms::Array{Float64} = [1.0, 2.0, 3.0, 4.0, 5.0]
	TMOSs = [["+", "+"], ["+", "-"], ["-", "-"], ["-", "+"]]


	confs = readdir(basename)
	println("confs: ", confs)
	# conf_new = findall(occursin.(r"^[0-9][0-9][0-9][0-9]_r[0-9]$", confs))
	conf_new = findall(occursin.("pion_mix_SIB", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	println("confs: ", confs)
	gammas = ["P5P5", "A1A1", "A2A2", "A3A3", "A4A4", "S0S0", "V1V1", "V2V2", "V3V3", "V4V4"]


	ncorr::Int32 = length(gammas) * (length(masses) * 2 * length(TMOSs) * (length(counterterms)))

	size::Int32 = ncorr * 2 * T #  ncorr *reim*T
	println("size: ", size)
	head = header(length(confs), T, L, ncorr, beta, kappa, masses, [+1.0, -1.0], [0.0], gammas, ["ll"], info_counterterms, counterterms, size)

	outfilename = "SIB_B48.dat"
	outfile = open(outfilename, "w")
	print(head, outfile)
	flush(outfile)

	println(head.Njack, " ", head.T)
	println(head.ncorr, " ", head.size)

	for (iconf, conf) in enumerate(confs)
		println(conf)
		write(outfile, Int32(iconf))
		hits = readdir(string(basename, "/", conf))

		# hits
		local conf_new = findall(occursin.(r"twop2_id[0-9]*_st[0-9]*\.h5", hits))
		hits_qcd = hits[conf_new]
		# hits_qcd = hits[conf_new[1:100]]
		# hits_qcd = hits[1:10]
		println("hits: ", length(hits_qcd))
		# println("hits: ", hits_qcd)

		### hits average 
		confname = string(basename, "/", conf, "/")
		corr = hits_SIB(outfile, confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)

		corr1::Array{Float64, 7} = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)

		# @time hits_average_SIB(outfile , corr1, confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)

		ig::Int = 1
		im = 1
		println("mass= ", masses[im], "   ", masses[im])
		println("GAMMA= ", gammas[ig])
		# println("corr", corr[im, 1, 1, ig, 1, 1:4, 1])
		# println("corr_S1", corr[im, 1, 1, ig, 2, 1:4, 1])
		# println("corr_S2", corr[im, 1, 1, ig, 3, 1:4, 1])
		# println("corr_P51", corr[im, 1, 1, ig, 4, 1:4, 1])
		# println("corr_P52", corr[im, 1, 1, ig, 5, 1:4, 1])
		ave = 0
		t = 1
		for i in 1:length(hits_qcd)
			ave += corr[i, im, 1, 1, ig, 1, t, 1]
		end
		ave /= length(hits_qcd)
		s = 0
		for i in 1:length(hits_qcd)
			s += (corr[i, im, 1, 1, ig, 1, t, 1] - ave)^2
		end
		s /= length(hits_qcd) - 1
		s = sqrt(s)
		println(ave / V, "  ", s / V)


		for ig in [1, 2, 6, 7]
			println("################################################################################################")
			println("gamma ", ig-1)
			iTMOS::Int = 1
			println("TM++")
			@printf("t     corr                       S_leg1_re                      S_leg2_re                     P5_leg1_im                P5_leg2_im\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V / sqrt(length(hits_qcd)))
			end
			@printf("t     corr_im                    S_leg1_im                     S_leg2_im                      P5_leg1_re               P5_leg2_re\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V / sqrt(length(hits_qcd)))
			end
			println("TM--")
			iTMOS = 3
			@printf("t     corr                       S_leg1_re                      S_leg2_re                     P5_leg1_im                P5_leg2_im\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V / sqrt(length(hits_qcd)))
			end
			@printf("t     corr_im                    S_leg1_im                     S_leg2_im                      P5_leg1_re               P5_leg2_re\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V / sqrt(length(hits_qcd)))
			end
			println("################################################################################################")
			iTMOS = 2
			println("OS+-")
			@printf("t     corr                       P5_leg1_re                      P5_leg2_re                     S_leg1_im                S_leg2_im\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V / sqrt(length(hits_qcd)))
			end
			@printf("t     corr_im                     P5_leg1_im                     P5_leg2_im                      S_leg1_re               S_leg2_re\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V / sqrt(length(hits_qcd)))
			end
			iTMOS = 4
			println("OS-+")
			@printf("t     corr                       P5_leg1_re                      P5_leg2_re                     S_leg1_im                S_leg2_im\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 2]) / V / sqrt(length(hits_qcd)))
			end
			@printf("t     corr_im                     P5_leg1_im                     P5_leg2_im                      S_leg1_re               S_leg2_re\n")
			for t in 1:T
				@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n", t - 1,
					mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 2, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V, std(corr[:, im, 1, iTMOS, ig, 3, t, 2]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 4, t, 1]) / V / sqrt(length(hits_qcd)),
					mean(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V, std(corr[:, im, 1, iTMOS, ig, 5, t, 1]) / V / sqrt(length(hits_qcd)))
			end
			println("################################################################################################")
		end

		conf_loop = findall(occursin.(r"twop3_id[0-9]*_st[0-9]*\.h5", hits))
		hits_loop = hits[conf_loop]
		# hits_qcd = hits[conf_new[1:100]]
		# hits_qcd = hits[1:10]
		println("hits: ", length(hits_qcd))
		counterterms_l::Array{Float64} = [1.0, 2.0, 3.0]
		TMOSs_l = [["+", "+"], ["-", "-"]]

		loop = ave_hits_loop(outfile, confname, head.T, hits_loop, masses, TMOSs_l, info_counterterms, counterterms_l, gammas)
		for ig in 1:16
			println("gamma ", ig-1)
			iTMOS::Int = 1
			@printf("bolla++:   %g    +I %g    \n", mean(loop[:, im, 1, iTMOS, ig, 1, 1])/V, mean(loop[:, im, 1, iTMOS, ig, 1, 2])/V )
			iTMOS = 2
			@printf("bolla--:   %g    +I %g    \n", mean(loop[:, im, 1, iTMOS, ig, 1, 1])/V, mean(loop[:, im, 1, iTMOS, ig, 1, 2])/V )

		end

		acc_data::Array{Float64, 5} = zeros(T, 16, 5, length(TMOSs), 2)
		read_open(conf, acc_data, basename, L, T, masses[1], TMOSs)

		for ig in 1:10
			println("################################################################################################")
			println("gamma ", ig-1)
			iTMOS=3
			println("TM ", TMOSs[iTMOS])
			for t in 1:T
				@printf("%-4d %-13.5g  %-13.5g      %-13.5g  %-13.5g\n", t-1, mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V,
				acc_data[t,ig,1,iTMOS,1],acc_data[t,ig,1,iTMOS,2]) 
			end
		end
		for ig in 1:10
			println("################################################################################################")
			println("gamma ", ig-1)
			iTMOS=4
			println("OS ", TMOSs[iTMOS])
			for t in 1:T
				@printf("%-4d %-13.5g  %-13.5g      %-13.5g  %-13.5g\n",t-1, mean(corr[:, im, 1, iTMOS, ig, 1, t, 1]) / V, mean(corr[:, im, 1, iTMOS, ig, 1, t, 2]) / V,
				acc_data[t,ig,1,iTMOS,1],acc_data[t,ig,1,iTMOS,2])  
			end
		end
		
		# for t in 1:T
		# 	@printf("%-4d %-13.5g    %-13.5g     %-13.5g  %-13.5g   %-13.5g\n", t - 1,
		# 		acc_data[t, ig,1, 1, 1],
		# 		acc_data[t, ig,2, 1, 1],
		# 		acc_data[t, ig,3, 1, 1],
		# 		acc_data[t, ig,4, 1, 2],
		# 		acc_data[t, ig,5, 1, 2])
		# end

		# @printf( "t     corr              S_leg1                    S_leg2                         P5_leg1                  P5_leg2\n" )
		# for t in 1:T
		# 	@printf("%-4d %-13.5g   %-13.5g     %-13.5g  %-13.5g   %-13.5g\n",t-1, 
		# 	corr1[im, 1, 1, ig, 1, t, 1],
		# 	corr1[im, 1, 1, ig, 2, t, 1],
		# 	corr1[im, 1, 1, ig, 3, t, 1],
		# 	corr1[im, 1, 1, ig, 4, t, 1],
		# 	corr1[im, 1, 1, ig, 5, t, 1] )
		# end

	end
	close(outfile)
end

main()
