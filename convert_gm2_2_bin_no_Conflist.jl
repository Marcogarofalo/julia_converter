using Printf

include("./read_hdf5.jl")
include("./gamma.jl")
# include("./read_LIBE_open.jl")
include("binning.jl")

include("modules_rew.jl")    # Load the file

using .modules_rew


function parse_data_file!(data::Array{Float64, 2}, filepath, Nconf::Int32, L::Int32)

    f = open(filepath, "r")
	current_header = ""
	for i in 1:Nconf
		line = readline(f)
        # println(line)
		if !isempty(line)
			error("non empty line where it should be")
		end
		line = readline(f)
        # println(line)
		if !startswith(line, "#")
			error("conf line does not start with #")
		end
		line = readline(f)
        # println(line)
		if !isempty(line)
			error("non empty line where it should be")
		end
		for t in 1:(L+1)
			line = readline(f)
			data[i, t] = parse(Float64, line)
		end
	end
    close(f)

end

function main()
	if length(ARGS)!=7
		println("usage: julia convert_sanfo_gmw.jl    gm2file binfile conflist   beta kappa  mu1 mu2")
		exit(1)
	end
	if ARGS[1]==ARGS[2]
		println("input and output must be different")
		exit(1)
	end

	println(ARGS[1])
	filename = ARGS[1]

	# Regex Explanation:
	# cB\.\d+\.([^_]+)  -> Matches cB. followed by any digits and another dot, then captures until '_' (96)
	# _r\..*_mu\.       -> Matches the constant _r., skips the middle text, and anchors at _mu.
	# ([^_]+)           -> Captures characters after _mu. until the next '_' (0.00072)
	# _([^_]+)\.txt     -> Captures everything between the last underscore and .txt (P5A0)
	pattern = r"c[a-zA-Z]\.\d+\.([^_]+)_r\..*_mu\.([^_]+)_([^_]+)\.txt"

	m = match(pattern, filename)
	if m !== nothing
		L = parse(Int32, m.captures[1])
		mu = parse(Float64, m.captures[2])
		gamma = m.captures[3]

		println("L: ", L)
		println("mu: ", mu)
		println("gamma: ", gamma)
	end
	num_lines = countlines(filename)
	println("Total lines: ", num_lines)
	if ((num_lines % (L+4))!=0)
		error("something is wrong num_lines ", num_lines, "   can not be divede by L+4 = ", L+4)
	end
	Nconfs::Int32 = div(num_lines, L+4)
	println(Nconfs)
	T::Int32 = L * 2

	data = zeros(Float64, Nconfs, T)

	parse_data_file!(data, filename, Nconfs, L)

    ### conf list
    # confs = readlines(ARGS[3])
	# if length(confs) != Nconfs
	# 	error("File length mismatch: from bin file $Nconfs confs, but conflist has  $(length(confs)) lines.")
	# end
    confs = fill("boh", Nconfs)
	

	T = L*2
	ncorr::Int32 = 1
	beta = parse(Float64, ARGS[4])
	kappa = parse(Float64, ARGS[5])
	masses::Vector{Float64} = [parse(Float64, ARGS[6]), parse(Float64, ARGS[7])]
	oranges::Vector{Float64} = []
	gamma_list::Vector{String} = [gamma]
	
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T
	mus::Array{Float64} = masses
    if(masses[1]!=mu || masses[2]!=mu)
        error("masses do not match")
    end

	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses, [0], [0.0], gamma_list, confs, [], oranges, sizeblock)

    conf_int::Vector{Int32} = Vector{Int32}(undef, length(confs))
	
	for (i, conf) in enumerate(confs)
        num_str = split(conf, '_')[1] 
		conf_int[i] = i
	end

    #### write 

	outfile = open(ARGS[2], "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)
	for (ic, conf) in enumerate(confs)
		write(outfile, conf_int[ic])
		for t in 1:T
			write(outfile, data[ic, t])
			write(outfile, Float64(0.0))
		end


	end
end

main()
