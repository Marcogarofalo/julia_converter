using HDF5
using Printf
import Base.write
using Profile
using DelimitedFiles

include("./read_hdf5.jl")
include("./gamma.jl")
include("./read_LIBE_open.jl")
include("binning.jl")


function char_before_match(s::String, pattern::String)
	pos = findfirst(pattern, s)
	pos1 = pos[1]
	mys = Vector{String}(undef, 2)  # Define a vector of strings with two elements
	if pos == nothing || pos[1] == 1
		error("Pattern not found or found at the beginning of the string")
	else
		rep::Char = s[pos1-1]
		mys[1] = replace(s[(pos1+1):length(s)], "gradflow.00" => "")
		if rep == 'a'
			mys[2] = mys[1] * "_r0"
		elseif rep == 'b'
			mys[2] = mys[1] * "_r1"
		elseif rep == 'c'
			mys[2] = mys[1] * "_r2"
		elseif rep == 'd'  # This is the last case
			mys[2] = mys[1] * "_r3"
		else
			error("Pattern not found or found at the beginning of the string")
		end

		return mys
	end
end

function main()

	if length(ARGS) != 2
		println("usage: julia convert_one_libe.jl  file   input.jl")
		exit(1)
	end
	outname::String = ARGS[1]
	include(ARGS[2])

	if length(masses) != length(name_quarks)
		println("error:")
		println("masses = ", masses)
		println("name_quarks = ", name_quarks)
		exit(1)
	end

	gamma_list::Vector{String} = ["gen_-Ig1g2",
		"gen_-Ig1g3",
		"gen_-Ig2g3",
		"gen_-Ig4g1",
		"gen_-Ig4g2",
		"gen_-Ig4g3",
		"gen_1",
		"gen_g1",
		"gen_g2",
		"gen_g3",
		"gen_g4",
		"gen_g5",
		"gen_g5g1",
		"gen_g5g2",
		"gen_g5g3",
		"gen_g5g4",
		"std_-Ig1g2",
		"std_-Ig1g3",
		"std_-Ig2g3",
		"std_-Ig4g1",
		"std_-Ig4g2",
		"std_-Ig4g3",
		"std_1",
		"std_g1",
		"std_g2",
		"std_g3",
		"std_g4",
		"std_g5",
		"std_g5g1",
		"std_g5g2",
		"std_g5g3",
		"std_g5g4",
	]

	## files
	files_n::Vector{String} = Vector{String}(undef, length(gamma_list) * length(name_quarks))
	count::Int32 = 1
	for j in 1:length(name_quarks)
		for i in 1:length(gamma_list)
			files_n[count] = basename_in * "/" * name_quarks[j] * "/new_" * gamma_list[i] * ".txt"
			count += 1
		end
	end
	fps::Vector{IOStream} = Vector{String}(undef, length(files_n))
	for (i, f) in enumerate(files_n)
		fps[i] = open(f)
	end
	println(files_n)


	## check file 1
	confs::Vector{String} = []
	for (i, line) in enumerate(eachline(fps[1]))
		if (occursin(Regex("^ # [0-9]*_r[0-9]*"), line))
			s = split(line, ' ')
			push!(confs, s[3])
		end
	end
	println(confs)
	println(length(confs))
	## reset
	close(fps[1])
	fps[1] = open(files_n[1])

	## set conf as integer
	conf_int::Vector{Int32} = Vector{Int32}(undef, length(confs))
	for (i, conf) in enumerate(confs)
		conf_int[i] = parse(Int32, replace(conf, Regex("_r[0-9]*") => ""))
	end

	# header
	ncorr::Int32 = length(gamma_list) * length(name_quarks)
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T

	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses, [0], [0.0], gamma_list, confs, [], [], sizeblock)


	outfile = open(outname, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)

	corr::Array{Float64, 3} = zeros(Float64, length(fps), T, 2)

	for (ic, conf) in enumerate(confs)

		for (fi, fp) in enumerate(fps)
			a = readline(fp)# empty
			a = readline(fp)# # 0000_r
			a = readline(fp)# empty
			# println("read head ", files_n[fi])
			for t in 1:T
				s = split(readline(fp), ' ')
				# corr[fi, t, :] = parse(readline(fp), Vector{Float64})
				corr[fi, t, 1] = parse(Float64, s[1])
				corr[fi, t, 2] = parse(Float64, s[2])
				# println(t-1,"  ",corr[fi, t, 1], corr[fi, t, 2], " ", files_n[fi])
			end
			# println("read empty")
			a = readline(fp)# empty 
		end

		write(outfile, conf_int[ic])
		for i in 1:length(fps)
			for t in 1:T
				write(outfile, corr[i, t, 1])
				write(outfile, corr[i, t, 2])
			end
		end


	end


end

main()
