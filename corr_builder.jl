using HDF5
using Printf
import Base.write
using Profile

function exact_exact(basename::String, T::Int32, nvec::Int32, confs::Array{String, 1})
	raw_data::Array{Float64, 3} = Array{Float64, 3}(undef, 2, nvec, nvec)
	for (iconf, conf) in enumerate(confs)
		filename = basename * "/" * conf * "/twop_exact_exact_nvec" * string(nvec) * ".h5"
		fid = h5open(filename, "r")
		println(keys(fid)[1])
		raw_data = read(fid, "littleD")
		# println(raw_data[1,1,1],"  ",raw_data[1,1,2])
		# println(raw_data[1,2,1],"  ",raw_data[1,2,2])
		# println(raw_data[2,1,1],"  ",raw_data[2,1,2])
		# println(raw_data[2,2,1],"  ",raw_data[2,2,2])
		# exit(1)
		close(fid)
	end
end


function stoch_stoch(basename::String, T::Int32, nvec::Int32, confs::Array{String, 1}, m::Float64, TMOSs::Vector{Vector{String}})
	# {128, 1, 4, 4, 4, 4, 2}
	raw_data::Array{Float64, 7} = Array{Float64, 7}(undef, T, 1, 4, 4, 4, 4, 2)
	for (iconf, conf) in enumerate(confs)


		hits = readdir(string(basename, "/", conf))
		pattern = string("twop_nev", nvec, "_id.._st...\\.h5")
		local conf_new = findall(occursin.(Regex(pattern), hits))
		hits = hits[conf_new]
		println(hits)
		for (i, hit) in enumerate(hits)
            filename = string(basename, "/", conf, "/", hit)
            fid = h5open(filename, "r")
            println(keys(fid)[1])
			for (iTMOS, TMOS) in enumerate(TMOSs)
				group::String = @sprintf("/%s/mesons/%c%.4e_%c%.4e_stoch_stoch", keys(fid)[1], TMOS[1], m, TMOS[2], m)
				# println(keys(fid)[1])
				raw_data = read(fid, group)
                # println(raw_data[1,1,1,1,1,1,1])
				# exit(1)
				# println(raw_data[1,1,1],"  ",raw_data[1,1,2])
				# println(raw_data[1,2,1],"  ",raw_data[1,2,2])
				# println(raw_data[2,1,1],"  ",raw_data[2,1,2])
				# println(raw_data[2,2,1],"  ",raw_data[2,2,2])
			end
            close(fid)
		end

	end

end

function main()
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B64/pion_defl/"
	T::Int32 = 128
	L::Int32 = 64
	beta::Float64 = 1.778000000000 ##check
	kappa::Float64 = 0.139426500000 ##check
	mass::Float64 = 7.2000e-04
	TMOSs = [["+", "+"], ["+", "-"]]
	nvec::Int32 = 400

	##find confs
	println(basename)
	confs = readdir(basename)
	conf_new = findall(occursin.("_r", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))
	# exact_exact(basename, T, nvec, confs)
	stoch_stoch(basename, T, nvec, confs, mass, TMOSs)

end

main()
