using HDF5
using Printf
import Base.write
using Profile

include("./read_hdf5.jl")
# qui i correlatori hanno 5 componenti
# le combinazioni sono 00, 01, 10, 02, 20
# dove 0 e' il propagator normale, 1 e' D^2 e 2 e' D g5 D

function main()
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_SIB"

	T::Int32 = 96
	L::Int32 = 48
	beta = 1.778000000000 ##check
	kappa = 0.139426500000 ##check
	masses = [1.8200e-02, 1.5000e-03, 3.0000e-03, 7.2000e-04]

	info_counterterms::Vector{Int32} = [5]
	# e::Float64 = 1.0000e-03 # before there was written by mistake 1e-2  
	counterterms::Array{Float64} = [1.0, 2.0, 3.0, 4.0, 5.0]
	TMOSs = [["+", "+"], ["+", "-"]]


	confs = readdir(basename)
	println("confs: ", confs)
	conf_new = findall(occursin.(r"^[0-9][0-9][0-9][0-9]_r[0-9]$", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	gammas = ["P5P5", "A1A1", "A2A2", "A3A3", "A4A4", "V1V1", "V2V2", "V3V3", "V4V4"]


	ncorr::Int32 = length(gammas) * (length(masses) * 2 * length(TMOSs) * (length(counterterms)))

	size::Int32 = ncorr * 2 * T #  ncorr *reim*T
	println("size: ", size)
	head = header(length(confs), T, L, ncorr, beta, kappa, masses, [+1.0, -1.0], [0.0], gammas, ["ll"], info_counterterms, counterterms, size)

	outfilename = "SIB_B48.dat"
	outfile = open(outfilename, "w")
	print(head, outfile)
	flush(outfile)

	println(head.Njack, " ", head.T)
	println(head.ncorr, " ", head.size)

	for (iconf, conf) in enumerate(confs)
		println(conf)
		write(outfile, Int32(iconf))
		hits = readdir(string(basename, "/", conf))

		# hits
		local conf_new = findall(occursin.(r"twop2_id[0-9]*_st[0-9]*\.h5", hits))
		hits_qcd = hits[conf_new]
		println("hits: ", length(hits_qcd))

		corr = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)
		### hits average 
		confname = string(basename, "/", conf, "/")
		@time hits_average_SIB(outfile, corr, confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)

		# write_hits_average(outfile, corr, confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)
		ig::Int = 1
		im = 1
		println("mass= ", masses[im], "   ", masses[im])
		println("GAMMA= ", gammas[ig])
		println("corr", corr[im, 1, 1, ig, 1, 1:4, 1])
		println("corr_S1", corr[im, 1, 1, ig, 2, 1:4, 1])
		println("corr_S2", corr[im, 1, 1, ig, 3, 1:4, 1])
		println("corr_P51", corr[im, 1, 1, ig, 4, 1:4, 2])
		println("corr_P52", corr[im, 1, 1, ig, 5, 1:4, 2])
		im = 4
        println("mass= ", masses[im], "   ", masses[im])
		println("GAMMA= ", gammas[ig])
		println("corr", corr[im, 1, 1, ig, 1, 1:4, 1])
		println("corr_S1", corr[im, 1, 1, ig, 2, 1:4, 1])
		println("corr_S2", corr[im, 1, 1, ig, 3, 1:4, 1])
		println("corr_P51", corr[im, 1, 1, ig, 4, 1:4, 1])
		println("corr_P52", corr[im, 1, 1, ig, 5, 1:4, 1])

		# println("00",corr[1, 1, 1, 1, 5, 1:4, 1])
		# println("ee -2 (00) + (-e-e)",corr[1, 1, 1, 1, 9, 1:4, 1]-2*corr[1, 1, 1, 1, 5, 1:4, 1]  + corr[1, 1, 1, 1, 1, 1:4, 1] )

		# println("-e0",corr[1, 1, 1, 1, 4, 1:4, 1])
		# println("e0",corr[1, 1, 1, 1, 6, 1:4, 1])
		# println("00",corr[1, 1, 1, 1, 5, 1:4, 1])
		# println("e0 -2 (00) + (-e0)",corr[1, 1, 1, 1, 6, 1:4, 1]-2*corr[1, 1, 1, 1, 5, 1:4, 1]  + corr[1, 1, 1, 1, 4, 1:4, 1] )

		# println("0-e",corr[1, 1, 1, 1, 2, 1:4, 1])
		# println("0e",corr[1, 1, 1, 1, 8, 1:4, 1])
		# println("00",corr[1, 1, 1, 1, 5, 1:4, 1])
		# println("0e -2 (00) + (0-e)",corr[1, 1, 1, 1, 8, 1:4, 1]-2*corr[1, 1, 1, 1, 5, 1:4, 1]  + corr[1, 1, 1, 1, 2, 1:4, 1] )
		# println("(ee -(-ee)- (e-e)+(-e,-e))/4", (corr[1, 1, 1, 1, 9, 1:4, 1]-corr[1, 1, 1, 1, 3, 1:4, 1] -corr[1, 1, 1, 1, 7, 1:4, 1] + corr[1, 1, 1, 1, 1, 1:4, 1] )/4)

		# println("OS: e0 -2 (00) + (-e0)",corr[1, 1, 2, 1, 6, 1:4, 1]-2*corr[1, 1, 2, 1, 5, 1:4, 1]  + corr[1, 1, 2, 1, 4, 1:4, 1] )

		# println("OS: 0e -2 (00) + (0-e)",corr[1, 1, 2, 1, 8, 1:4, 1]-2*corr[1, 1, 2, 1, 5, 1:4, 1]  + corr[1, 1, 2, 1, 2, 1:4, 1] )

	end
	close(outfile)
end

main()
