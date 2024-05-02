using HDF5
using Printf
import Base.write
using Profile

include("./read_hdf5.jl")
include("./gamma.jl")
include("./read_LIBE_open.jl")
include("binning.jl")


function binning(corr_all::Array{Float64, 8}, head::header, Nconfs::Int64, TMOSs::Vector{Vector{String}})
	masses = head.mus
	gammas = head.gammas
	counterterms = head.oranges
	T = head.T
	Nb::Int = head.Njack
	m1list::Array{Float64, 1} = [masses[1], masses[1]]
	slice::Array{Float64, 1} = Array{Float64, 1}(undef, Nconfs)
	corr_bin::Array{Float64, 8} = zeros(Float64, Nb, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)
	bin_slice::Array{Float64, 1} = Array{Float64, 1}(undef, Nb)
	for (im, m) in enumerate(masses)
		m1list[2] = m
		for (im1, m1) in enumerate(m1list)
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for ig in 1:(length(gammas))
					for ic in 1:(length(counterterms))
						for t in 1:T
							for reim in 1:2
								slice .= corr_all[:, im, im1, iTMOS, ig, ic, t, reim]
								bin_slice = bin_intoN(slice, Nb)
								corr_bin[:, im, im1, iTMOS, ig, ic, t, reim] .= bin_slice
							end
						end
					end
				end
			end
		end
	end
	return (corr_bin)
end

function main()
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE"
	T::Int32 = 96
	L::Int32 = 48
	beta = 1.778000000000 ##check
	kappa = 0.139426500000 ##check
	masses = [1.8250e-02, 6.8400e-03, 5.0400e-03, 3.6000e-03, 2.1600e-03, 7.2000e-04]
	Nb = 50

	info_counterterms::Vector{Int32} = [9, 5]
	e::Float64 = 1.0000e-03 # before there was written by mistake 1e-2  
	counterterms = [-e, 0, e, -e, 0, e, -e, 0, e, 0, 1, 2, 3, 4]
	TMOSs = [["+", "+"], ["+", "-"]]

	if (sum(info_counterterms)!=length(counterterms))
		prinln("Error: sum(info_counterterms): ",sum(info_counterterms), "  differs from length(counterterms): ",length(counterterms))
		exit(2)
	end
	confs = readdir(basename)

	conf_new = findall(occursin.(r"[0-9][0-9][0-9][0-9]_r", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	gamma_list::Vector{String} = ["G5",#1
		"Id",#2
		"G5x",#3
		"G5y",#4
		"G5z",#5
		"G5t",#6
		"Gx",#7
		"Gy",#8
		"Gz",#9
		"Gt",#10
		"sxy", "sxz", "syz", "stx", "sty", "stz"]
	###### open
	# gammas::Vector{String} = Vector{String}(undef, length(gamma_list) * length(gamma_list))
	# for (i, gi) in enumerate(gamma_list)
	# 	for (j, gj) in enumerate(gamma_list)
	# 		gammas[(i)+(j-1)*length(gamma_list)] = gi * gj
	# 	end
	# end
	###### else
	gammas::Vector{String} = gamma_list

	# println(gammas)

	g5GI_list::Vector{gamma_struct} = [G5 * G5,#1
		G5 * Id,#2
		G5 * G5 * Gx,#3
		G5 * G5 * Gy,#4
		G5 * G5 * Gz,#5
		G5 * G5 * Gt,#6
		G5 * Gx,#7
		G5 * Gy,#8
		G5 * Gz,#9
		G5 * Gt,#10
		G5 * sxy, G5 * sxz, G5 * syz, G5 * stx, G5 * sty, G5 * stz]

	# for (i,g) in enumerate(g5GI_list)
	# 	g5GI_list[i] = g
	# end

	ncorr::Int32 = length(gammas) * (length(masses) * 2 * length(TMOSs) * (length(counterterms)))

	size::Int32 = ncorr * 2 * T #  ncorr *reim*T
	println("size: ", size)
	println("Nconfs: ", length(confs))
	# consider only configurations with data inside
	confs_with_data::Vector{Int}  = []
	for (iconf, conf) in enumerate(confs)
		hits::Vector{String} = readdir(string(basename, "/", conf))
		pattern::String = string("^twop_id[0-9]*_st[0-9]*\\.h5\$")
		local conf_new = findall(occursin.(Regex(pattern), hits))
		hits = hits[conf_new]
		if ( length(hits)> 0 )
			push!(confs_with_data, iconf)
		end
	end
	confs = confs[confs_with_data]
	println("Nconfs with data: ", length(confs))

	head::header = header(Nb, T, L, ncorr, beta, kappa, masses, [+1.0, -1.0], [0.0], gammas, ["ll"], info_counterterms, counterterms, size)

	outfilename = "LIBE_B48.dat"
	outfile = open(outfilename, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)

	println("Njack  T = ", head.Njack, " ", head.T)
	println("Ncorr size = ", head.ncorr, " ", head.size)
	corr_all::Array{Float64, 8} = zeros(Float64, length(confs), length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)

	# for (iconf, conf) in enumerate(confs)
	# Threads.@threads
	for iconf in 1:length(confs)
		conf = confs[iconf]
		# write(outfile, Int32(iconf))
		hits = readdir(string(basename, "/", conf))

		# hits

		# println("hits: ", hits_qcd)

		##### open
		# corr::Array{Float64, 7} = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)
		# fill!(corr, 0.0)
		# @time read_LIBE_open(conf, corr, basename, L, T, masses, TMOSs, info_counterterms, g5GI_list)
		# corr_all[iconf, :, :, :, :, :, :, :] .= corr[:, :, :, :, :, :, :]
		### else
		@time read_LIBE!(conf, corr_all, iconf, basename, L, T, masses, TMOSs, info_counterterms, g5GI_list)
		flush(stdout)

	end


	### binning
	@time corr_bin = binning(corr_all, head, length(confs), TMOSs)

	# writing 
	@time m1list::Array{Float64, 1} = [masses[1], masses[1]]
	for ni in 1:Nb
		write(outfile, Int32(ni))
		for (im, m) in enumerate(masses)
			m1list[2] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)
					for ig in 1:(length(gammas))
						for ic in 1:(length(counterterms))
							for t in 1:T
								for reim in 1:2
									write(outfile, corr_bin[ni, im, im1, iTMOS, ig, ic, t, reim])
								end
							end
						end
					end
				end
			end
		end
	end
	close(outfile)
end

main()
