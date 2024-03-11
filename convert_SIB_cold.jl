using HDF5
using Printf
import Base.write
using Profile
using Statistics

include("./read_hdf5.jl")
# qui i correlatori hanno 5 componenti
# le combinazioni sono 00, 01, 10, 02, 20
# dove 0 e' il propagator normale, 1 e' D^2 e 2 e' D g5 D

function hits_SIB(outfile::IOStream, 
	confname::String, T::Int32, hits_qcd::Vector{String},
	masses::Array{Float64, 1}, TMOSs::Vector{Vector{String}},
	info_counterterms::Array{Int32, 1}, counterterms::Array{Float64, 1},
	gammas::Vector{String})

	# h5ls -r file gives: {96, 1, 3, 3, 16, 2}
	raw_data::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, 5, 1, T)
	m1list::Array{Float64, 1}=[masses[1],masses[1]]
	corr::Array{Float64,8} = zeros(Float64, length(hits_qcd), length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)

	for (ihit, hit) in enumerate(hits_qcd)
		fname::String = string(confname, hit)
		fid = h5open(fname, "r")
		for (im, m) in enumerate(masses)
			m1list[1]=m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)

					combo::String = @sprintf("%s/mesons/%c%.4e_%c%.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					# raw_data = fid[combo][:, id, 1, :]
					raw_data = read(fid, combo)
					# raw_data = h5read(fname, combo)

					for (ig, g) in enumerate(gammas)
						id::Int32 = id_of_gamma[ig]
						for ic in 1:length(counterterms)
							
							for t in 1:T
								corr[ihit, im, im1, iTMOS, ig, ic, t, 1] += raw_data[1, id, ic, 1, t]
								corr[ihit, im, im1, iTMOS, ig, ic, t, 2] += raw_data[2, id, ic, 1, t]
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
	return(corr)
end
function main()
	# basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_SIB"
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/cold"

	T::Int32 = 8
	L::Int32 = 4
	V::Int = L*L*L
	beta = 1.778000000000 ##check
	kappa = 0.14 ##check
	masses = [0.001]

	info_counterterms::Vector{Int32} = [5]
	# e::Float64 = 1.0000e-03 # before there was written by mistake 1e-2  
	counterterms::Array{Float64} = [1.0, 2.0, 3.0, 4.0, 5.0]
	TMOSs = [["+", "+"], ["+", "-"]]


	confs = readdir(basename)
	println("confs: ", confs)
	# conf_new = findall(occursin.(r"^[0-9][0-9][0-9][0-9]_r[0-9]$", confs))
	conf_new = findall(occursin.("pion_mix_SIB", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	println("confs: ", confs)
	gammas = ["P5P5", "A1A1", "A2A2", "A3A3", "A4A4", "V1V1", "V2V2", "V3V3", "V4V4"]


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
		println("hits: ", hits_qcd)

		### hits average 
		confname = string(basename, "/", conf, "/")
		corr=hits_SIB(outfile , confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)

		corr1::Array{Float64,7} = zeros(Float64,  length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)

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
		ave=0
		t=1
		for i in 1:length(hits_qcd)
			ave +=corr[i,im, 1, 1, ig, 1, t, 1]
		end
		ave/=length(hits_qcd)
		s=0
		for i in 1:length(hits_qcd)
			s +=(corr[i,im, 1, 1, ig, 1, t, 1]-ave)^2
		end
		s/=length(hits_qcd)-1
		s=sqrt(s)
		println(ave/V, "  ",s/V);
		@printf( "t     corr                       S_leg1_re                      S_leg2_re                      P5_leg1_im                P5_leg2_im\n" )
		for t in 1:T
			@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n",t-1, 
			mean(corr[:,im, 1, 1, ig, 1, t, 1])/V,std(corr[:,im, 1, 1, ig, 1, t, 1])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 2, t, 1])/V,std(corr[:,im, 1, 1, ig, 2, t, 1])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 3, t, 1])/V,std(corr[:,im, 1, 1, ig, 3, t, 1])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 4, t, 2])/V,std(corr[:,im, 1, 1, ig, 4, t, 2])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 5, t, 2])/V,std(corr[:,im, 1, 1, ig, 5, t, 2])/V/sqrt(length(hits_qcd)) )
		end
		@printf( "t     corr_im                     S_leg1_im                     S_leg2_im                      P5_leg1_re               P5_leg2_re\n" )
		for t in 1:T
			@printf("%-4d %-13.5g( %-10.5g)    %-13.5g( %-10.5g)     %-13.5g( %-10.5g)  %-13.5g( %-10.5g)   %-13.5g( %-10.5g)\n",t-1, 
			mean(corr[:,im, 1, 1, ig, 1, t, 2])/V,std(corr[:,im, 1, 1, ig, 1, t, 2])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 2, t, 2])/V,std(corr[:,im, 1, 1, ig, 2, t, 2])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 3, t, 2])/V,std(corr[:,im, 1, 1, ig, 3, t, 2])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 4, t, 1])/V,std(corr[:,im, 1, 1, ig, 4, t, 1])/V/sqrt(length(hits_qcd)),
			mean(corr[:,im, 1, 1, ig, 5, t, 1])/V,std(corr[:,im, 1, 1, ig, 5, t, 1])/V/sqrt(length(hits_qcd)))
		end

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
