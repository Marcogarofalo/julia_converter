using HDF5
using Printf
import Base.write
using Profile

include("./read_hdf5.jl")
include("./gamma.jl")
include("./read_LIBE_open.jl")
include("binning.jl")

include("modules_rew.jl")    # Load the file

using .modules_rew     


function load_input(in)
	include(in)
	# pattern_after_rep = isdefined(Main, :pattern_after_rep) ? Main.pattern_after_rep : "/"
	# len_rep_name = isdefined(Main, :len_rep_name) ? Main.len_rep_name : 0
	# rep_a_in_number = isdefined(Main, :rep_a_in_number) ? Main.rep_a_in_number : 0
	
	# return basename_in,  T, L, beta, kappa, masses, outname, pattern_after_rep, rep_a_in_number, len_rep_name,  confs
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
		println("usage: julia script.jl     input.jl")
		exit(1)
	end
	
	# outname::String = ARGS[1]
	basename_in,  T, L, beta, kappa, masses, outname, pattern_after_rep, rep_a_in_number, len_rep_name,  confs = load_input(ARGS[1])

	gamma_list::Vector{String} = ["t", "P", "Eplaq", "Esym", "tsqEplaq", "tsqEsym", "Wsym", "Qsym"]
	println("pattern_after_rep: ", pattern_after_rep	)
	ncorr::Int32 = length(gamma_list)
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T

	corr::Array{Float64, 2} = zeros(Float64, 8, T)

	conf_formated::Vector{String} = Vector{String}(undef, length(confs))
	conf_int::Vector{Int32} = Vector{Int32}(undef, length(confs))
	two_s = Vector{String}(undef, 2)
	for (i, conf) in enumerate(confs)
		two_s = char_before_match(conf, pattern_after_rep, rep_a_in_number, len_rep_name)
		conf_int[i] = parse(Int32, two_s[1])
		conf_formated[i] = two_s[2]
	end

	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses, [0], [0.0], gamma_list, conf_formated, [], [], sizeblock)

	outfile = open(outname, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)

	errors::Int = 0
	for (ic, conf) in enumerate(confs)
		filename = basename_in * "/" * conf
		if !isfile(filename)
			println("Error: ", filename, "  not found")
			println( conf_formated[ic])
			errors += 1
		end
	end
	if errors != 0
		println("confs missing : ", errors, "  ")
		exit(1)
	end

	for (ic, conf) in enumerate(confs)

		filename = basename_in * "/" * conf
		fp = open(filename)
		a = readline(fp)
		println(conf)
		Nt::Int32 = 0
		for (i, line) in enumerate(eachline(fp))
			(a, b, c, d, e, f, g, e, f) = NTuple{9}(eachsplit(line, ' '))
			id = parse(Int32, a)
			corr[1, i] = parse(Float64, b)
			corr[2, i] = parse(Float64, c)
			corr[3, i] = parse(Float64, d)
			corr[4, i] = parse(Float64, e)
			corr[5, i] = parse(Float64, f)
			corr[6, i] = parse(Float64, g)
			corr[7, i] = parse(Float64, e)
			corr[8, i] = parse(Float64, f)
			Nt += 1
		end
		if (Nt != T)
			println("Error: Nt != T , Nt = ", Nt, " T = ", T)
			exit(1)
		end
		write(outfile, conf_int[ic])
		for i in 1:8
			for t in 1:T
				write(outfile, corr[i, t])
				write(outfile, Float64(0.0))
			end
		end


	end
	println(length(confs), " configurations processed.")
	println("Output written to ", outname)

end

main()
