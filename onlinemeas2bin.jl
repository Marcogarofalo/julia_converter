using HDF5
using Printf
import Base.write
using Profile
using DelimitedFiles
using Statistics: Statistics
using Profile

include("./read_hdf5.jl")
include("./gamma.jl")
include("./read_LIBE_open.jl")
include("binning.jl")

include("modules_rew.jl")    # Load the file

using .modules_rew

function read_file!(corr, filename_hit::String, head::header)
	# println("filename_hit: ", filename_hit)
	f = open(filename_hit, "r")
	line_count::Int32 = 0
	for lines in readlines(f)
		line_count += 1
	end
	if (line_count != div((head.T + 2) * head.ncorr, 2))
		println("Error: ", filename_hit, "  line counted:", line_count, " != T*Ncorr/2: ", div((head.T + 2) * head.ncorr, 2))
		exit(1)
	end
	# close(f)
	# f = open(filename_hit, "r")
	seekstart(f) # = seek(f, 0)

	for icor in 1:head.ncorr

		line = readline(f)
		s = split(line, r"\s+")
		corr[icor, 1] += parse(Float64, s[4])
		for t in 2:(div(head.T, 2))
			line = readline(f)
			s = split(line, r"\s+")
			corr[icor, t] += parse(Float64, s[4])
			corr[icor, head.T-t+2] += parse(Float64, s[5])
		end
		line = readline(f)
		s = split(line, r"\s+")
		corr[icor, div(head.T + 2, 2)] += parse(Float64, s[4])

		# for t in 1:(head.T)
		# 	println("t: ", t-1, "  corr[", icor, ", ", t, "] = ", corr[icor, t])
		# end
	end
	close(f)

end

function load_input(in)
	include(in)
	# len_rep_name = isdefined(Main, :len_rep_name) ? Main.len_rep_name : 0
	# rep_a_in_number = isdefined(Main, :rep_a_in_number) ? Main.rep_a_in_number : 0
	# Nhits = isdefined(Main, :Nhits) ? Main.Nhits : 24
	# return basename_in, monomial, T, L, beta, kappa, masses_in, masses_out, outname, pattern_after_rep, rep_a_in_number, len_rep_name, ave_sources, confs, Nhits

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

	Nhits = Base.invokelatest(isdefined, @__MODULE__, :Nhits) ? Base.invokelatest(getfield, @__MODULE__, :Nhits) : 24	
	ave_sources = Base.invokelatest(isdefined, @__MODULE__, :ave_sources) ? Base.invokelatest(getfield, @__MODULE__, :ave_sources) : error("ave_sources not defined in input file")	


	return basename_in, monomial, T, L, beta, kappa, masses_in, masses_out, outname, pattern_after_rep, rep_a_in_number, len_rep_name, ave_sources, confs, Nhits

end

function main()

	if length(ARGS) != 1
		println("usage: julia convert_one_libe.jl   input.jl")
		exit(1)
	end

	basename_in, monomial, T, L, beta, kappa, masses_in, masses_out, outname, pattern_after_rep, rep_a_in_number, len_rep_name, ave_sources, confs, Nhits = load_input(ARGS[1])


	gamma_list::Vector{String} = ["P5P5", "A0P5", "A0A0"]

	## files
	file_n::Array{String} = Vector{String}(undef, length(confs))
	count::Int32 = 1
	for j in 1:length(confs)
		file_n[count] = basename_in * "/" * confs[j]
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
	replica_dir::Vector{String} = Vector{String}(undef, length(confs))
	for (i, conf) in enumerate(confs)
		conf_int[i] = parse(Int32, replace(conf, Regex(".*/onlinemeas\\.s...\\.") => ""))
		confs_name[i] = char_before_match(conf, pattern_after_rep, rep_a_in_number, len_rep_name)[2]
		replica_dir[i] = replace(conf, Regex("/onlinemeas\\.s.*") => "")
		# println(conf_int[i], " ", confs_name[i], " ", replica_dir[i])
	end


	# header
	ncorr::Int32 = length(gamma_list)
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T

	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses_in, [0], [0.0], gamma_list, confs_name, [], masses_out, sizeblock)


	outfile = open(outname, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)


	wU::Vector{Float64} = Vector{Float64}(undef, length(confs))
	corr::Array{Float64, 2} = zeros(Float64, ncorr, T)


	errors::Int = 0
	for (ic, conf) in enumerate(confs)

		if (ave_sources)
			conf_in_dir::Vector{String} = readdir(string(basename_in, "/", replica_dir[ic]))
			conf06 = @sprintf("%06d", conf_int[ic])
			pattern::String = string("^onlinemeas\\.s...\\." * conf06 * "\$")
			# println(pattern)
			hits = findall(occursin.(Regex(pattern), conf_in_dir))
			if length(hits) != Nhits
				println("Error: ", conf, "  no Nhits hits found, it has ", length(hits))
				errors += 1
			end
		end
	end

	if errors != 0
		println("confs missing : ", errors, "  ")
		exit(1)
	end

	for (ic, conf) in enumerate(confs)

		if (ave_sources)


			# regex_name=replace(conf, Regex("onlinemeas.s....") => "onlinemeas.s....")
			# println(regex_name) 
			conf_in_dir::Vector{String} = readdir(string(basename_in, "/", replica_dir[ic]))
			conf06 = @sprintf("%06d", conf_int[ic])
			pattern::String = string("^onlinemeas\\.s...\\." * conf06 * "\$")
			# println(pattern)
			hits = findall(occursin.(Regex(pattern), conf_in_dir))
			println("conf: ", conf, "  Nhits: ", length(hits))
			if length(hits) != Nhits
				println("Error: ", conf, "  no Nhits hits found")
				exit(1)
			end
			for hit in hits
				filename_hit = string(basename_in, "/", replica_dir[ic], "/", conf_in_dir[hit])
				# println("filename_hit: ", filename_hit)
				read_file!(corr, filename_hit, head)
			end

			for icorr in 1:ncorr
				for t in 1:T
					corr[icorr, t] /= length(hits)
				end
			end

			# print(conf_in_dir[hits])
			# exit(1)
			# pattern::String = string("^twop_id[0-9]*_st[0-9]*\\.h5\$")
			# local conf_new = findall(occursin.(Regex(pattern), hits))
		else
			println("do something when we do not need to average the sources: ")
		end

		write(outfile, conf_int[ic])
		for icorr in 1:ncorr
			for t in 1:T
				write(outfile, corr[icorr, t])
				write(outfile, 0)
				corr[icorr, t] = 0.0
			end
		end

		# f = open(files_n[ic], "r")
		# line_count::Int32 = 0
		# for lines in readlines(f)
		# 	# s = split(lines, ' ')
		# 	# corr[ic, t, 1] = parse(Float64, s[7])
		# 	# increment line_count
		# 	line_count += 1

		# 	# print the line
		# 	# println(lines)        
		# end
		# close(f)
		# f = open(files_n[ic], "r")
		# ws::Vector{Float64} = Vector{Float64}(undef, line_count)
		# for (i, lines) in enumerate(readlines(f))
		# 	s = split(lines, ' ')
		# 	ws[i] = parse(Float64, s[7])
		# 	# println(ws[i])
		# end

		# # ws = ws[1:(div(line_count,4))]
		# wU[ic] = compute_wU(ws, monomial)

		# close(f)
	end

	# ave::Float64 = Statistics.mean(wU)
	# corr::Array{Float64, 3} = zeros(Float64, length(wU), T, 2)
	# # normalization
	# for (ic, conf) in enumerate(wU)
	# 	# corr[ic, 1, 1] = length(wU) / (sum(exp.(wU .- wU[ic])))
	# 	corr[ic, 1, 1] = exp(wU[ic] - ave )
	# end

	# ## effective configurations
	# bw = compute_wU(wU, TM_monomial());
	# bw2 = compute_wU(wU .*2, TM_monomial());
	# Nfactor = 1.0 / exp(bw2-2*bw)
	# Nfactor = Nfactor / length(wU);
	# println("Nfactor: ", Nfactor)

	# # sorted_indices = sortperm(corr[:, 1, 1])
	# # println(sorted_indices)
	# # println(corr[sorted_indices, 1, 1])
	# for (ic, conf) in enumerate(wU)
	# 	# println(ic, "  ",wU[ic], "  ",corr[ic, 1, 1], "    ", sum(corr[:, 1, 1]))
	# 	write(outfile, conf_int[ic])
	# 	for t in 1:T
	# 		write(outfile, wU[ic])
	# 		write(outfile, corr[ic, t, 1])
	# 		# println("W(U) = ",wU[ic], "  exp(W(U)-<W(U)>) = ", corr[ic, t, 1])
	# 	end
	# end
	# # println(confs[694])
	# # println(confs[635])
	# # println(confs[103])

end


main() # Compile once

# # Run for profiling
# Profile.clear()
# @profile main()

# # Write the results to a file named "profile_output.txt"
# open("profile_output.txt", "w") do io
#     Profile.print(io, format=:flat, sortedby=:count) # You can customize the options here
# end

# println("Profile results written to profile_output.txt")