using HDF5
using Printf
import Base.write
using Profile
using DelimitedFiles
using Statistics: Statistics

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
		mys[1] = replace(s[(pos1+length(pattern)):length(s)], "reweighting_factor.data." => "")
		mys[1] = @sprintf("%04d", parse(Int32, mys[1]))
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

struct TM_monomial end
struct OS_monomial end
monomial_type = Union{TM_monomial, OS_monomial}

function compute_wU(ws::Vector{Float64}, monomial::TM_monomial)
	# println("TM")
	m = Statistics.mean(ws)
	wU = m + log(sum(exp.(ws .- m)))
	return wU
end
function compute_wU(ws::Vector{Float64}, monomial::OS_monomial)
	# println("OS  ",length(ws),"  ws = "	,ws)
	m::Float64 = Statistics.mean(ws)
	# wU = m / 2.0 + log(sum(exp.((ws .- m) ./ 2.0))) # wrong formula
	wU::Float64 = (m + log(sum(exp.(ws .- m)))) / 2.0
	return wU
end

function main()

	if length(ARGS) != 1
		println("usage: julia convert_one_libe.jl   input.jl")
		exit(1)
	end
	# outname::String = ARGS[1]
	include(ARGS[1])

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
		confs_name[i] = char_before_match(conf, pattern_after_rep)[2]
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

	for (ic, conf) in enumerate(confs)
		f = open(files_n[ic], "r")
		line_count::Int32 = 0
		for lines in readlines(f)
			# s = split(lines, ' ')
			# corr[ic, t, 1] = parse(Float64, s[7])
			# increment line_count
			line_count += 1

			# print the line
			# println(lines)        
		end
		close(f)
		f = open(files_n[ic], "r")
		ws::Vector{Float64} = Vector{Float64}(undef, line_count)
		for (i, lines) in enumerate(readlines(f))
			s = split(lines, ' ')
			ws[i] = parse(Float64, s[7])
			# println(ws[i])
		end

		ws_red = ws[1:(div(line_count,reduce_sources_by))]
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
		println("conf: ", confs_name[ic], " ic-1: ", ic - 1, "  rU*e^{-ave}: ", corr[ic, 1, 1], "  wU: ", wU[ic], "  ave: ", ave)
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

end

main()
