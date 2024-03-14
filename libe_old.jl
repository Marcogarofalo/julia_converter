using HDF5
using Printf
import Base.write
using Profile

include("./read_hdf5.jl")

function main()
basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE_1"

	confs = readdir(basename)

	conf_new = findall(occursin.("_r", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))

	T::Int32 = 96
	L::Int32 = 48
	beta = 1.778000000000 ##check
	kappa = 0.139426500000 ##check

	masses = [1.8200e-02, 1.5000e-03, 3.0000e-03, 7.2000e-04]
	info_counterterms::Vector{Int32} = [4, 4, 4] # this tell ous the differences in the array of ds, 1= #dmu, 1=#dkappa, 2=#de
	counterterms = [+7.0000e-01, +7.0000e-02, -7.0000e-01, -7.0000e-02,
		+1.0000e-04, +1.0000e-05, -1.0000e-04, -1.0000e-05,
		1.0000e-02, 1.0000e-03, -1.0000e-02, -1.0000e-03]
	TMOSs = [["+", "+"], ["+", "-"]]

	gammas = ["P5P5", "A1A1", "A2A2", "A3A3", "A4A4", "V1V1", "V2V2", "V3V3", "V4V4"]
	# (length(counterterms) + 1) : +1 because there is the correlator without counterterms
	# length(masses) * 2 * length(TMOSs) : 2 because the masses combinations of the twopoint function: (m1,m1) and (mn,m1)
	ncorr::Int32 = length(gammas) * (length(masses) * 2 * length(TMOSs) * (length(counterterms) + 1))

	size::Int32 = ncorr * 2 * T #  ncorr *reim*T
	println(size)
	head = header(length(confs), T, L, ncorr, beta, kappa, masses, [+1.0, -1.0], [0.0], gammas, ["ll"], info_counterterms, counterterms, size)
	outfilename = "LIBE_B48.dat"
	outfile = open(outfilename, "w")
	print(head, outfile)

	flush(outfile)
	in = open(outfilename, "r")
	read(in, head.Njack)
	read(in, head.T)
	println(head.Njack, " ", head.T)


	# confs = confs[[1, 2]]
	for (iconf, conf) in enumerate(["1680_r1"])#enumerate(confs)
		println(conf)
		write(outfile, Int32(iconf))
		hits = readdir(string(basename, "/", conf))

		##only ending with .h5
		# println(hits)
		local conf_new = findall(occursin.(r".*?\.h5$", hits))
		hits = hits[conf_new]
		# println(hits)

		# local conf_new = findall(.!occursin.("old", hits))
		# hits = hits[conf_new]
		# println(hits)


		# QCD hits
		local conf_new = findall(occursin.(r"twop_id.._st...\.h5", hits))
		hits_qcd = hits[conf_new]
		println("QCD hits: ", length(hits_qcd))

		corr = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), sum(info_counterterms) + 1, T, 2)
		# corr::Array{Float64,7}(0, length(masses), 2, length(TMOSs), length(gammas), sum(info_counterterms)+1 , T, 2)
		# corr::Array{Float64,7}(zero, length(masses), 2, length(TMOSs), length(gammas), sum(info_counterterms)+1 , T, 2)
		### hits average QCD
		confname = string(basename, "/", conf, "/")
		@time qcd_part(outfile, corr, confname, head.T, hits_qcd, masses, TMOSs, info_counterterms, counterterms, gammas)


		# qed_hits
		local conf_new = findall(occursin.(r"twop_id.._st..._qed....h5", hits))
		hits_qed = hits[conf_new]
		println("qed hits: ", length(hits_qed))

		# @time QED_part(outfile, corr, confname, head.T, hits_qed, masses, TMOSs, info_counterterms, counterterms, gammas)

		### print hits average
		@time write_hits_average(outfile, corr, confname, head.T, hits_qed, masses, TMOSs, info_counterterms, counterterms, gammas)
        
        println("ee: ",corr[1, 1, 1, 1, 9+1, 1:4, 1])
        # println("ee: ",corr[1, 1, 1, 1, 10+1, 1:4, 1])
        println("-e-e: ",corr[1, 1, 1, 1, 11+1, 1:4, 1])
        # println("ee: ",corr[1, 1, 1, 1, 12+1, 1:4, 1])
		println("00: ",corr[1, 1, 1, 1, 1, 1:4, 1])
		println("ee -2 (00) + (-e-e)",corr[1, 1, 1, 1, 9+1, 1:4, 1]-2*corr[1, 1, 1, 1, 1, 1:4, 1]  + corr[1, 1, 1, 1, 11+1, 1:4, 1] )
		println("e'e': ",corr[1, 1, 1, 1, 10+1, 1:4, 1])
        println("-e'-e': ",corr[1, 1, 1, 1, 12+1, 1:4, 1])
		println("e'e' -2 (00) + (-e'-e')",corr[1, 1, 1, 1, 10+1, 1:4, 1]-2*corr[1, 1, 1, 1, 1, 1:4, 1]  + corr[1, 1, 1, 1, 12+1, 1:4, 1] )

		println("strange ",masses[1])
		println("masses")
		println("dmu=",counterterms[2])
		println("dmudmu: ",corr[1, 1, 1, 1, 1+2, 1:4, 1])
		println("-dmu-dmu: ",corr[1, 1, 1, 1, 1+4, 1:4, 1])
		println("[dmudmu - (-dmu-dmu)]/2eps: ",(corr[1, 1, 1, 1, 1+2, 1:4, 1]-corr[1, 1, 1, 1, 1+4, 1:4, 1])/(2*counterterms[2]*masses[1])   )
		println("dmu=",counterterms[1])
		println("dmudmu: ",corr[1, 1, 1, 1, 1+1, 1:4, 1])
		println("-dmu-dmu: ",corr[1, 1, 1, 1, 1+3, 1:4, 1])
		println("[dmudmu - (-dmu-dmu)]/2eps: ",(corr[1, 1, 1, 1, 1+1, 1:4, 1]-corr[1, 1, 1, 1, 1+3, 1:4, 1])/(2*counterterms[1]*masses[1])   )


		println("kappa")
		println("dk=",counterterms[2+4])
		println("dkdk: ",corr[1, 1, 1, 1, 1+4+2, 1:4, 1])
		println("-dk-dk: ",corr[1, 1, 1, 1, 1+4+4, 1:4, 1])
		println("[dkdk - (-dk-dk)]/2eps: ",(corr[1, 1, 1, 1, 1+4+2, 1:4, 1]-corr[1, 1, 1, 1, 1+4+4, 1:4, 1])/(2*counterterms[4+2])   )
		println("dk=",counterterms[1+4])
		println("dkdk: ",corr[1, 1, 1, 1, 1+4+1, 1:4, 1])
		println("-dk-dk: ",corr[1, 1, 1, 1, 1+4+3, 1:4, 1])
		println("[dkdk - (-dk-dk)]/2eps: ",(corr[1, 1, 1, 1, 1+4+1, 1:4, 1]-corr[1, 1, 1, 1, 1+4+3, 1:4, 1])/(2*counterterms[4+1])   )

		println("light ",masses[4])
		println("masses")
		println("dmu=",counterterms[2])
		println("dmudmu: ",corr[4, 1, 1, 1, 1+2, 1:4, 1])
		println("-dmu-dmu: ",corr[4, 1, 1, 1, 1+4, 1:4, 1])
		println("[dmudmu - (-dmu-dmu)]/2eps: ",(corr[4, 1, 1, 1, 1+2, 1:4, 1]-corr[4, 1, 1, 1, 1+4, 1:4, 1])/(2*counterterms[2]*masses[4])   )
		println("dmu=",counterterms[1])
		println("dmudmu: ",corr[4, 1, 1, 1, 1+1, 1:4, 1])
		println("-dmu-dmu: ",corr[4, 1, 1, 1, 1+3, 1:4, 1])
		println("[dmudmu - (-dmu-dmu)]/2eps: ",(corr[4, 1, 1, 1, 1+1, 1:4, 1]-corr[4, 1, 1, 1, 1+3, 1:4, 1])/(2*counterterms[1]*masses[1])   )
		println("kappa")
		println("dk=",counterterms[2+4])
		println("dkdk: ",corr[4, 1, 1, 1, 1+4+2, 1:4, 1])
		println("-dk-dk: ",corr[4, 1, 1, 1, 1+4+4, 1:4, 1])
		println("[dkdk - (-dk-dk)]/2eps: ",(corr[4, 1, 1, 1, 1+4+2, 1:4, 1]-corr[4, 1, 1, 1, 1+4+4, 1:4, 1])/(2*counterterms[4+2])   )
		println("dk=",counterterms[1+4])
		println("dkdk: ",corr[4, 1, 1, 1, 1+4+1, 1:4, 1])
		println("-dk-dk: ",corr[4, 1, 1, 1, 1+4+3, 1:4, 1])
		println("[dkdk - (-dk-dk)]/2eps: ",(corr[4, 1, 1, 1, 1+4+1, 1:4, 1]-corr[4, 1, 1, 1, 1+4+3, 1:4, 1])/(2*counterterms[4+1])   )



	end
	close(outfile)
	# fid = h5open("/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE_1/1680_r1/twop_id00_st000_qed000.h5", "r")
	# pipi = read_correlator(fid, "/sx00sy00sz00st00/mesons/+1.5000e-03_+1.5000e-03_de+1.0000e-02_de+1.0000e-02", 1)
	# println(length(pipi))
end
#@profile main()

main()
