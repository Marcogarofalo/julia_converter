#using HDF5
using Printf
import Base.write
using Profile
using DelimitedFiles
using Statistics: Statistics

include("./read_hdf5.jl")
include("./gamma.jl")
#include("./read_LIBE_open.jl")
include("binning.jl")

include("modules_rew.jl")    # Load the file

using .modules_rew     


function load_input(in)
	println("Loading input file: ", in)
	# evalfile(in)
	include(in)
	# println(Base.invokelatest(getfield, @__MODULE__, :basename_in))
	# len_rep_name = isdefined(Main, :len_rep_name) ? Main.len_rep_name : 0
	# rep_a_in_number = isdefined(Main, :rep_a_in_number) ? Main.rep_a_in_number : 0
	# check_mult = isdefined(Main, :check_mult) ? Main.check_mult : true
	# reduce_sources_by = isdefined(Main, :reduce_sources_by) ? Main.reduce_sources_by : Int32(1)

	basename_in = Base.invokelatest(isdefined, @__MODULE__, :basename_in) ? Base.invokelatest(getfield, @__MODULE__, :basename_in) : error("basename_in not defined in input file")
	monomial = Base.invokelatest(isdefined, @__MODULE__, :monomial) ? Base.invokelatest(getfield, @__MODULE__, :monomial) : error("monomial not defined in input file")
	T = Base.invokelatest(isdefined, @__MODULE__, :T) ? Base.invokelatest(getfield, @__MODULE__, :T) : error("T not defined in input file")
	L = Base.invokelatest(isdefined, @__MODULE__, :L) ? Base.invokelatest(getfield, @__MODULE__, :L) : error("L not defined in input file")
	beta = Base.invokelatest(isdefined, @__MODULE__, :beta) ? Base.invokelatest(getfield, @__MODULE__, :beta) : error("beta not defined in input file")
	kappa = Base.invokelatest(isdefined, @__MODULE__, :kappa) ? Base.invokelatest(getfield, @__MODULE__, :kappa) : error("kappa not defined in input file")
	masses_in = Base.invokelatest(isdefined, @__MODULE__, :masses_in) ? Base.invokelatest(getfield, @__MODULE__, :masses_in) : error("masses_in not defined in input file")
	masses_out = Base.invokelatest(isdefined, @__MODULE__, :masses_out) ? Base.invokelatest(getfield, @__MODULE__, :masses_out) : error("masses_out not defined in input file")
	outname = Base.invokelatest(isdefined, @__MODULE__, :outname) ? Base.invokelatest(getfield, @__MODULE__, :outname) : error("outname not defined in input file")
	pattern_after_rep = Base.invokelatest(isdefined, @__MODULE__, :pattern_after_rep) ? Base.invokelatest(getfield, @__MODULE__, :pattern_after_rep) : error("pattern_after_rep not defined in input file")

	rep_a_in_number = Base.invokelatest(isdefined, @__MODULE__, :rep_a_in_number) ? Base.invokelatest(getfield, @__MODULE__, :rep_a_in_number) : 0	
	len_rep_name = Base.invokelatest(isdefined, @__MODULE__, :len_rep_name) ? Base.invokelatest(getfield, @__MODULE__, :len_rep_name) : 0	
	confs = Base.invokelatest(isdefined, @__MODULE__, :confs) ? Base.invokelatest(getfield, @__MODULE__, :confs) : error("confs not defined in input file")	



	return basename_in, monomial, T, L, beta, kappa, masses_in, masses_out, outname, pattern_after_rep, rep_a_in_number, len_rep_name,  confs
end


function main()

	if length(ARGS) != 1
		println("usage: julia convert_one_libe.jl   input.jl")
		exit(1)
	end
	
	
	# outname::String = ARGS[1]
	# include(ARGS[1])
	basename_in, monomial, T, L, beta, kappa, masses_in, masses_out, outname, pattern_after_rep, rep_a_in_number, len_rep_name,  confs = load_input(ARGS[1])

	gamma_list::Vector{String} = [""]

	## files
	files_n::Vector{String} = Vector{String}(undef, length(confs))
	count::Int32 = 1
	for j in 1:length(confs)
		files_n[count] = basename_in * "/" * confs[j]
		count += 1
	end
	# fps::Vector{IOStream} = Vector{String}(undef, length(files_n))
	# for (i, f) in enumerate(files_n)
	# 	fps[i] = open(f)
	# end
	# println(files_n)



	## set conf as integer
	conf_int::Vector{Int32} = Vector{Int32}(undef, length(confs))
	confs_name::Vector{String} = Vector{String}(undef, length(confs))
	for (i, conf) in enumerate(confs)
		conf_int[i] = parse(Int32, replace(conf, Regex(".*/reweighting_factor.data.") => ""))
		confs_name[i] = char_before_match(conf, pattern_after_rep, rep_a_in_number, len_rep_name)[2]
	end


	# header
	ncorr::Int32 = length(gamma_list)
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T

	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses_in, [0], [0.0], gamma_list, confs_name, [], masses_out, sizeblock)


	outfile = open(outname, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)

	errors::Int = 0
	for (ic, conf) in enumerate(confs)
		if !isfile(files_n[ic])
				println("Error: ",files_n[ic] , "  not found")
				errors += 1
		end
	end
	if errors != 0
		println("confs missing : ", errors, " from a list of : ", length(confs))
		exit(1)
	end

	wU::Vector{Float64} = Vector{Float64}(undef, length(confs))
    line_count::Vector{Int32} = zeros(Int32, length(confs))
	for (ic, conf) in enumerate(confs)
		f = open(files_n[ic], "r")

		for lines in readlines(f)
			# s = split(lines, ' ')
			# corr[ic, t, 1] = parse(Float64, s[7])
			# increment line_count
			line_count[ic] += 1

			# print the line
			# println(lines)        
		end
		close(f)
		f = open(files_n[ic], "r")
		ws::Vector{Float64} = Vector{Float64}(undef, line_count[ic])
		for (i, lines) in enumerate(readlines(f))
			s = split(lines, ' ')
			ws[i] = parse(Float64, s[7])
			kappa = parse(Float64, s[3])
			if (abs(kappa - head.kappa) > 1e-10)
				error("kappa in input ", head.kappa, " not equal of the one of file ", files_n[ic], " : ", kappa)
			end
			if (check_mult)
				mu1 = parse(Float64, s[6]) / (2.0 * kappa)
				if (abs(mu1 - head.mus[1]) > 1e-10)
					error("mu in input denominator ", head.mus[1], " not equal of the one of file ", files_n[ic], " : ", mu1)
				end
				kappa = parse(Float64, s[4])
				mu2 = parse(Float64, s[5]) / (2.0 * kappa)
				if (abs(mu2 - head.oranges[1]) > 1e-10)
					error("mu in input numerator ", head.oranges[1], " not equal of the one of file ", files_n[ic], " : ", mu2)
				end
			end
			# println(ws[i])
		end

		ws_red::Vector{Float64} = -ws[1:(div(line_count[ic], reduce_sources_by))]
		wU[ic] = compute_wU(ws_red, monomial)
		# println("conf: ", confs_name[ic], " ic-1: ", ic - 1,  "  wU:",wU[ic])

		close(f)
	end
	# deleteat!(wU, 694)
	# deleteat!(wU, 635)
	# deleteat!(wU, 103)
	# println(exp.(wU))
	ave::Float64 = Statistics.mean(wU)
	corr::Array{Float64, 3} = zeros(Float64, length(wU), T, 2)
	# normalization
	for (ic, conf) in enumerate(wU)
		# corr[ic, 1, 1] = length(wU) / (sum(exp.(wU .- wU[ic])))
		corr[ic, 1, 1] = exp(wU[ic] - ave)
		println("conf: ", confs_name[ic], " hits: ", line_count[ic], " hits used: ", div(line_count[ic], reduce_sources_by), "  rU*e^{-ave}: ", corr[ic, 1, 1], "  wU: ", wU[ic], "  ave: ", ave)
	end

	## effective configurations
	bw = compute_wU(wU, TM_monomial())
	bw2 = compute_wU(wU .* 2, TM_monomial())
	Nfactor = 1.0 / exp(bw2 - 2 * bw)
	Nfactor = Nfactor / length(wU)
	println("Nfactor: ", Nfactor)

	# sorted_indices = sortperm(corr[:, 1, 1])
	# println(sorted_indices)
	# println(corr[sorted_indices, 1, 1])
	for (ic, conf) in enumerate(wU)
		# println(ic, "  ",wU[ic], "  ",corr[ic, 1, 1], "    ", sum(corr[:, 1, 1]))
		write(outfile, conf_int[ic])
		for t in 1:T
			write(outfile, wU[ic])
			write(outfile, corr[ic, t, 1])
			# println("W(U) = ",wU[ic], "  exp(W(U)-<W(U)>) = ", corr[ic, t, 1])
		end
	end
	# println(confs[694])
	# println(confs[635])
	# println(confs[103])
	close(outfile)
	println(length(wU), " configurations processed.")
	println("Output written to ", outname)
end

main()
