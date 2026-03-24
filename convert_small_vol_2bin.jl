#print(map(x->string(x, x), ARGS))
using Printf

include("./read_hdf5.jl")
include("./gamma.jl")
# include("./read_LIBE_open.jl")
include("binning.jl")

include("modules_rew.jl")    # Load the file

using .modules_rew

function main()
	println(length(ARGS))
	if length(ARGS)!=7
		println("usage: julia convert_sanfo_gmw.jl  binfile  gm2file conflist   beta kappa  mu1 mu2")
		exit(1)
	end
	if ARGS[1]==ARGS[2]
		println("input and output must be different")
		exit(1)
	end

	println("reading: ", ARGS[1])
	fs=filesize(ARGS[1])
	println("size: ", fs)
	io = open(ARGS[1], "r");

	Nconfs = read(io, Int64)
	Nhits = read(io, Int64)
	T::Int32 = read(io, Int64)
	Nsub::Int32 = read(io, Int64)

	mu1 = read(io, Float64)
	mu2 = read(io, Float64)


	println("Nconfs: ", Nconfs)
	println("Nhits: ", Nhits)
	println("T: ", T)
	println("Nsub: ", Nsub)
	confs = fill("boh", Nconfs)
	for i in 1:Nconfs
		# Read 7 bytes at once and convert to a String
		# Note: Julia uses 1-based indexing
		v = String(read(io, 7))

		# Append to the collection
		confs[i] = v
		

		# Print the element followed by spaces
		# print(confs[i], "  ")
	end
	# println() # Final newline

	# if (fs-sizeof(Int64)*4-(Nconfs*Nsub*sizeof(Float64)*(T/2+1)) != 0)
	# 	println("size ", fs, " and Nconfs ", Nconfs, "x", Nsub, "  do not match")
	# 	exit(1)
	# end

	chunk_size=div(T,2)+1
	#data = zeros(Float64, (Nconfs,chunk_size ))
	# data=Array{Float32}(undef, chunk_size, Nconfs)
	data::Array{Float64} = zeros(Float64, Nconfs, T)

	for i in 1:Nconfs
		for t in 1:chunk_size
			for k in 1:Nsub
				data[i, t]+=read(io, Float64)
			end
			data[i, t]/=Float64(Nsub)
			# println(data[i,t])
		end
	end
	# wrong order to read it in one go
	# read!(io, data)


	L = div(T,2)
	ncorr::Int32 = 1
	beta = parse(Float64, ARGS[4])
	kappa = parse(Float64, ARGS[5])
	masses::Vector{Float64} = [mu1, mu2]
	oranges::Vector{Float64} = []
	gamma_list::Vector{String} = ["P5P5"]
	close(io)
	sizeblock::Int32 = ncorr * 2 * T #  ncorr *reim*T
	mus::Array{Float64} = masses
	head::header = header(length(confs), T, L, ncorr, beta, kappa, masses, [0], [0.0], gamma_list, confs, [], oranges, sizeblock)


	conf_int::Vector{Int32} = Vector{Int32}(undef, length(confs))

	for (i, conf) in enumerate(confs)
		num_str = split(conf, '_')[1]
		conf_int[i] = parse(Int32, num_str)
	end
	#### write 
	println("writing ",ARGS[2])
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
	# out = open(ARGS[2], "w");
	# for i in 1:Nconfs
	#     println(out,"")
	#     @printf(out,"# %04d_r0\n",i)
	#     println(out,"")
	#     for t in 1:chunk_size
	#         println(out,data[i,t])
	#         # println(out,data[i,t])
	#     end
	# end

end

main()

