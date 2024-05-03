using HDF5
using Printf
import Base.write
using Profile

include("./read_hdf5.jl")
include("binning.jl")


function binning(corr_all::Array{Float64, 8}, head::header, Nconfs::Int64, TMOSs::Vector{Float64})
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
	if length(ARGS) != 1
		println("usage: julia convert_one_libe.jl  Nb")
		exit(1)
	end
	basename = "B48"
	Nb::Int32 = parse(Int32, ARGS[1])
	confs::Vector{String} = readdir(string(basename, "/"))
	conf_new = findall(occursin.(r"[0-9][0-9][0-9][0-9]_r", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	outfilename = "LIBE_B48.dat"


    ## init 
	conf = confs[1]
    file = open("B48/" * conf, "r")
    head::header = read_header(file)
    corr_all::Array{Float64, 8} = zeros(Float64, length(confs), length(head.mus), 2, length(head.rs), length(head.gammas), length(head.oranges), head.T, 2)
    close(file)

	for iconf in 1:length(confs)
		conf = confs[iconf]

		file = open("B48/" * conf, "r")
		local head::header = read_header(file)

		n::Int32 = read(file, Int32)

		##### open
		m1list::Array{Float64, 1} = [head.mus[1], head.mus[1]]
		for (im, m) in enumerate(head.mus)
			m1list[2] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(head.rs)
					for ig in 1:(length(head.gammas))
						for ic in 1:(length(head.oranges))
							for t in 1:head.T
								for reim in 1:2
									corr_all[iconf, im, im1, iTMOS, ig, ic, t, reim]=read(file, Float64)
								end
							end
						end
					end
				end
			end
		end
		close(file)

	end


	headw::header = header(Nb, head.T, head.L, head.ncorr, head.beta, head.kappa, head.mus, head.rs, head.thetas, head.gammas, head.smearing, head.bananas, head.oranges, head.size)
	### binning
	@time corr_bin = binning(corr_all, head, length(confs), head.rs)


	outfile = open(outfilename, "w")
	print(head, outfile)
	flush(outfile)
	flush(stdout)
	# writing 
	m1list = [head.mus[1], head.mus[1]]
	for ni in 1:Nb
		write(outfile, Int32(ni))
		for (im, m) in enumerate(head.mus)
			m1list[2] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(head.rs)
					for ig in 1:(length(head.gammas))
						for ic in 1:(length(head.oranges))
							for t in 1:head.T
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
