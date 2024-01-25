using HDF5
using Printf
import Base.write
using Profile




function read_correlator(raw_data, fid, dataset, id)
	id = id

	# buf = Vector{Float64}(undef,length(fid[dataset]))
	# h5d_read(fid[dataset], memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf)
	# raw_data=fid[dataset][:,id,1,:]

	# buf=read(fid,dataset=>Float64)
	# raw_data=@view  read(fid,dataset=>Float64)[:,id,1,:]




end

struct header
	Njack::Int32
	T::Int32
	L::Int32
	ncorr::Int32
	beta::Float64
	kappa::Float64
	mus::Array{Float64}
	rs::Array{Float64}
	thetas::Array{Float64}
	gammas::Array{String}
	smearing::Array{String}
	bananas::Array{Int32}
	oranges::Array{Float64}
	size::Int32
end

function get_st(in)
	start_string = "st"
	end_string = "_qed"
	s = string("(?<=", start_string, ").*?(?=", end_string, ")")
	pattern = Regex(s)
	match_result = match(pattern, in)
	if match_result !== nothing
		return match_result.match
	else
		s = string("(?<=", start_string, ").*?(?=_)")
		pattern = Regex(s)
		match_result = match(pattern, in)
	end
	return match_result.match
end

function mywrite(file, mus::Array{Float64})
	n::Int32 = length(mus)
	write(file, htol(n))
	write(file, htol(mus))
end

function mywrite(file, mus::Array{Int32})
	n::Int32 = length(mus)
	write(file, htol(n))
	write(file, htol(mus))
end

function mywrite(file, gammas::Array{String})
	n::Int32 = length(gammas)
	write(file, htol(n))
	for g in gammas
		for char in g
			write(file, char)
		end
		write(file, '\0')
	end
end


function print(head::header, file)

	write(file, htol(head.Njack))
	write(file, htol(head.T))
	write(file, htol(head.L))
	write(file, htol(head.ncorr))
	write(file, htol(head.beta))
	write(file, htol(head.kappa))

	mywrite(file, head.mus)
	mywrite(file, head.rs)
	mywrite(file, head.thetas)
	mywrite(file, head.gammas)
	mywrite(file, head.smearing)
	mywrite(file, head.bananas)
	mywrite(file, head.oranges)
	write(file, head.size)

end

const basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE_1"

#"pseudoscalar, scalar, g5g1, g5g2, g5g3, g5g4, g1, g2, g3, g4,s12,s13,s23,s41,s42,s43"
const id_of_gamma::Array{Int32, 1} = [1, 3, 4, 5, 6, 7, 8, 9, 10]
#["P5P5","A1A1","A2A2","A3A3","A4A4","V1V1","V2V2","V3V3","V4V4"]


function qcd_part(outfile, corr, confname::String, T::Int32, hits_qcd::Vector{String}, masses::Array{Float64, 1}, TMOSs::Vector{Vector{String}},
	info_counterterms::Array{Int32, 1}, counterterms::Array{Float64, 1}, gammas::Vector{String})
	# raw_data = Array{Float64, 2}(undef, 2, T)
	raw_data = Array{Float64, 4}(undef, 2, 16, 1, T)
	# corr = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), 1 + info_counterterms[1] + info_counterterms[2], T, 2) #= combination same mass, first mass =#
	for (ihit, hit) in enumerate(hits_qcd)
		fname = string(confname, hit)
		fid = h5open(fname, "r")
		for (im, m) in enumerate(masses)
			for (im1, m1) in enumerate([m, masses[1]])
				for (iTMOS, TMOS) in enumerate(TMOSs)

					combo::String = @sprintf("%s/mesons/%c%.4e_%c%.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					# raw_data = fid[combo][:, id, 1, :]
					raw_data = read(fid, combo)
					# memorymap
					# raw_data = fid[combo]
					# raw_data = HDF5.readmmap(raw_data)

					for (ig, g) in enumerate(gammas)
						id = id_of_gamma[ig]
						for t in 1:T
							corr[im, im1, iTMOS, ig, 1, t, 1] += raw_data[1, id, 1, t]
							corr[im, im1, iTMOS, ig, 1, t, 2] += raw_data[2, id, 1, t]
						end
					end
					## dmu 
					for i in 1:info_counterterms[1]
						offsave = 1 + i
						combo = @sprintf("%s/mesons/%c%.4e_%c%.4e_dmu%+-.4e_dk%+-.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1, counterterms[i], 0.0)
						# raw_data = fid[combo][:, id, 1, :]
						raw_data = read(fid, combo)
						# raw_data = fid[combo]
						# raw_data = HDF5.readmmap(raw_data)

						for (ig, g) in enumerate(gammas)
							id = id_of_gamma[ig]
							for t in 1:T
								corr[im, im1, iTMOS, ig, offsave, t, 1] += raw_data[1, id, 1, t]
								corr[im, im1, iTMOS, ig, offsave, t, 2] += raw_data[2, id, 1, t]
							end
						end
					end
					## dk
					for i in 1:info_counterterms[2]
						ipo = i + info_counterterms[1]
						offsave = 1 + info_counterterms[1] + i
						combo = @sprintf("%s/mesons/%c%.4e_%c%.4e_dmu%+-.4e_dk%+-.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1, 0.0, counterterms[ipo])
						# raw_data = fid[combo][:, id, 1, :]
						raw_data = read(fid, combo)
						# raw_data = fid[combo]
						# raw_data = HDF5.readmmap(raw_data)

						for (ig, g) in enumerate(gammas)
							id = id_of_gamma[ig]
							for t in 1:T
								corr[im, im1, iTMOS, ig, offsave, t, 1] += raw_data[1, id, 1, t]
								corr[im, im1, iTMOS, ig, offsave, t, 2] += raw_data[2, id, 1, t]
							end
						end
					end


				end
			end
		end
		close(fid)
	end

	# corr_qcd ./= length(hits_qcd)

	# writing
	for (im, m) in enumerate(masses)
		for (im1, m1) in enumerate([m, masses[1]])
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for (ig, g) in enumerate(gammas)

					for i in 1:(1+info_counterterms[1]+info_counterterms[2])
						for t in 1:T
							corr[im, im1, iTMOS, ig, i, t, 1] /= length(hits_qcd)
							corr[im, im1, iTMOS, ig, i, t, 2] /= length(hits_qcd)
							# write(outfile, corr_qcd[im, im1, iTMOS, ig, i, t, 1])
							# write(outfile, corr_qcd[im, im1, iTMOS, ig, i, t, 2])

						end
					end
				end
			end
		end
	end

end


function QED_part(outfile, corr, confname::String, T::Int32, hits_qed::Vector{String}, masses::Array{Float64, 1}, TMOSs::Vector{Vector{String}},
	info_counterterms::Array{Int32, 1}, counterterms::Array{Float64, 1}, gammas::Vector{String})
	# corr = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), info_counterterms[3], T, 2)
	raw_data = Array{Float64, 2}(undef, 2, T)
	for (ihit, hit) in enumerate(hits_qed)
		# println(hit)
		fname = string(confname, hit)

		fid = h5open(fname, "r")
		combo = "combo string"
		for (im, m) in enumerate(masses)
			for (im1, m1) in enumerate([m, masses[1]])
				for (iTMOS, TMOS) in enumerate(TMOSs)
					for (ig, g) in enumerate(gammas)

						## dk
						for i in 1:info_counterterms[3]
							ipo = i + info_counterterms[1] + info_counterterms[2]
							offsave = 1 + ipo
							combo = @sprintf("%s/mesons/%c%.4e_%c%.4e_de%+-.4e_de%+-.4e", keys(fid)[1], TMOS[1], m, TMOS[2], m1, counterterms[ipo], counterterms[ipo])

							id = id_of_gamma[ig]

							raw_data = fid[combo][:, id, 1, :]
							for t in 1:T
								corr[im, im1, iTMOS, ig, offsave, t, 1] += raw_data[1, t]
								corr[im, im1, iTMOS, ig, offsave, t, 2] += raw_data[2, t]
							end
						end
					end
				end
			end
		end

		close(fid)
	end


	# corr ./= length(hits_qed)
	### end hits average

	for (im, m) in enumerate(masses)
		for (im1, m1) in enumerate([m, masses[1]])
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for (ig, g) in enumerate(gammas)

					for i in (1+info_counterterms[1]+info_counterterms[2]):info_counterterms[3]
						for t in 1:T
							corr[im, im1, iTMOS, ig, i, t, 1] /= length(hits_qed)
							corr[im, im1, iTMOS, ig, i, t, 2] /= length(hits_qed)
							# write(outfile, corr[im, im1, iTMOS, ig, i, t, 1])
							# write(outfile, corr[im, im1, iTMOS, ig, i, t, 2])
						end
					end
				end

			end
		end
	end
end


function write_hits_average(outfile, corr, confname, T, hits_qed, masses, TMOSs, info_counterterms, counterterms, gammas)
	for (im, m) in enumerate(masses)
		for (im1, m1) in enumerate([m, masses[1]])
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for (ig, g) in enumerate(gammas)

					for i in 1:(sum(info_counterterms)+1)
						for t in 1:T
							# corr[im, im1, iTMOS, ig, i, t, 1] /=hits_qed
							# corr[im, im1, iTMOS, ig, i, t, 2] /=hits_qed
							write(outfile, corr[im, im1, iTMOS, ig, i, t, 1])
							write(outfile, corr[im, im1, iTMOS, ig, i, t, 2])
						end
					end
				end

			end
		end
	end
end

function main()

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
	for (iconf, conf) in enumerate(confs)
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

		@time QED_part(outfile, corr, confname, head.T, hits_qed, masses, TMOSs, info_counterterms, counterterms, gammas)

		### print hits average
		@time write_hits_average(outfile, corr, confname, head.T, hits_qed, masses, TMOSs, info_counterterms, counterterms, gammas)
		

	end
	close(outfile)
	# fid = h5open("/leonardo_scratch/large/userexternal/sbacchio/B48/pion_mix_LIBE_1/1680_r1/twop_id00_st000_qed000.h5", "r")
	# pipi = read_correlator(fid, "/sx00sy00sz00st00/mesons/+1.5000e-03_+1.5000e-03_de+1.0000e-02_de+1.0000e-02", 1)
	# println(length(pipi))
end
#@profile main()

main()
# end
## final id in c
# counterterm = none+ dmu1 + dmu2 + ... + dk1+dk2+...  + de1 + de2 + ..... 
# id = counterterm + Ncounterter *(TMOS + 2 *( m1 + 2*(m ) ))
# m1= [m, m[0]]