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

	if length(ARGS) != 2
		println("usage: julia convert_one_libe.jl  file   input.jl")
		exit(1)
	end

	
	include(ARGS[2])

	info_counterterms::Vector{Int32} = [9, 5]
	counterterms = [-e, 0, e, -e, 0, e, -e, 0, e, 0, 1, 2, 3, 4]

	if (sum(info_counterterms) != length(counterterms))
		prinln("Error: sum(info_counterterms): ", sum(info_counterterms), "  differs from length(counterterms): ", length(counterterms))
		exit(2)
	end

	conf = ARGS[1]
	println("confs: ", conf)
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
		# "sxy", "sxz", "syz", "stx", "sty", "stz"
	]
	#### open diag, P5*any,any*p5
	gammas::Vector{String} = Vector{String}(undef, length(gamma_list) * 3)
	for (i, gi) in enumerate(gamma_list)
		gammas[i] = gi * gi
	end
	for (i, gi) in enumerate(gamma_list)
		gammas[i+length(gamma_list)] = gi * gamma_list[1]
	end
	for (i, gi) in enumerate(gamma_list)
		gammas[i+2*length(gamma_list)] = gamma_list[1]* gi
	end
	###### open
	# gammas::Vector{String} = Vector{String}(undef, length(gamma_list) * length(gamma_list))
	# for (i, gi) in enumerate(gamma_list)
	# 	for (j, gj) in enumerate(gamma_list)
	# 		gammas[(i)+(j-1)*length(gamma_list)] = gi * gj
	# 	end
	# end
	###### else
	# gammas::Vector{String} = gamma_list

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
		G5 * Gt, #10
		# G5 * sxy, G5 * sxz, G5 * syz, G5 * stx, G5 * sty, G5 * stz
	]

	g5GI_source_sink::Array{gamma_struct, 2} = Array{gamma_struct, 2}(undef, 3 * length(g5GI_list), 2)
	println(size(g5GI_source_sink))
	for (i, gi) in enumerate(g5GI_list)
		g5GI_source_sink[i, 1] = gi
		g5GI_source_sink[i, 2] = gi
	end
	for (i, gi) in enumerate(g5GI_list)
		g5GI_source_sink[i+length(g5GI_list), 1] = gi
		g5GI_source_sink[i+length(g5GI_list), 2] = g5GI_list[1]
	end
	for (i, gi) in enumerate(g5GI_list)
		g5GI_source_sink[i+2*length(g5GI_list), 1] = g5GI_list[1]
		g5GI_source_sink[i+2*length(g5GI_list), 2] = gi
	end

	if (size(g5GI_source_sink)[1] != length(gammas))
		println("Error gamma arrray are not the same")
	end

	# for (i,g) in enumerate(g5GI_list)
	# 	g5GI_list[i] = g
	# end

	ncorr::Int32 = length(gammas) * (length(masses) * 2 * length(TMOSs) * (length(counterterms)))

	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T
	println("sizeblock: ", sizeblock)
	println("Nconfs: ", length(conf))
	# consider only configurations with data inside

	hits::Vector{String} = readdir(string(basename_in, "/", conf))
	pattern::String = string("^twop_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new]
	if (length(hits) == 0)
		println("non hits for this conf: ", conf)
		exit(1)
	end


	head::header = header(Nb, T, L, ncorr, beta, kappa, masses, [+1.0, -1.0], [0.0], gammas, ["ll"], info_counterterms, counterterms, sizeblock)


	outfile = open(  outname * "/" * conf, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)

	println("Njack  T = ", head.Njack, " ", head.T)
	println("Ncorr sizeblock = ", head.ncorr, " ", head.size)


	##### open
	corr::Array{Float64, 7} = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)
	@time read_LIBE_open(conf, corr, basename_in, L, T, masses, TMOSs, info_counterterms, g5GI_source_sink)

	# writing 
	@time m1list::Array{Float64, 1} = [masses[1], masses[1]]
	write(outfile, Int32(1))
	for (im, m) in enumerate(masses)
		m1list[2] = m
		for (im1, m1) in enumerate(m1list)
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for ig in 1:(length(gammas))
					for ic in 1:(length(counterterms))
						for t in 1:T
							for reim in 1:2
								write(outfile, corr[ im, im1, iTMOS, ig, ic, t, reim])
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
