using HDF5
using Printf
import Base.write
using Profile
using FFTW

include("./gamma.jl")

function symm_t!(data::Array{Float64, 4})
	T::Int = size(data)[1]

	for t in 2:div(T, 2)
		for iTMOS in 1:2
			for i in 1:10
				for reim in 1:2
					data[t, i, iTMOS, reim] += data[(T-t+2), i, iTMOS, reim]
					data[t, i, iTMOS, reim] /= 2.0
					data[(T-t+2), i, iTMOS, reim] = data[t, i, iTMOS, reim]
				end
			end
		end
	end
end


function Gamma_contraction(acc_data::Array{Float64, 4}, raw_data::Array{Float64, 7}, T::Int, iTMOS::Int)
	for t in 1:T
		# for i in 1:10
		# 	for a in 1:4
		# 		for b in 1:4
		# 			for c in 1:4
		# 				for d in 1:4

		# 					tmp::Complex{Float64} = (raw_data[1, a, b, c, d, 1, t] + (raw_data[2, a, b, c, d, 1, t])im) * gamma[i][a, d] * gamma[i][c, b]
		# 					acc_data[t, i, iTMOS, 1] += tmp.re
		# 					acc_data[t, i, iTMOS, 2] += tmp.im
		# 				end
		# 			end
		# 		end
		# 	end
		# end
		for i in 1:10
			for j in 1:4
				for k in 1:4
					tmp::Complex{Float64} = (raw_data[1, j, Gamma[i].col[k], k, Gamma[i].col[j], 1, t] + (raw_data[2, j, Gamma[i].col[k], k, Gamma[i].col[j], 1, t])im) * Gamma[i].val[j] * Gamma[i].val[k]
					acc_data[t, i, iTMOS, 1] += tmp.re
					acc_data[t, i, iTMOS, 2] += tmp.im
				end
			end
		end

	end
end

function stoch_stoch(conf::String, acc_data::Array{Float64, 4}, basename::String, L::Int, T::Int, nvec::Int, m::Float64, TMOSs::Vector{Vector{String}})
	# {128, 1, 4, 4, 4, 4, 2}
	raw_data::Array{Float64, 7} = Array{Float64, 7}(undef, 2, 4, 4, 4, 4, 1, T)


	hits = readdir(string(basename, "/", conf))
	pattern = string("^twop_nev", nvec, "_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new]
	println("Nhits:   ", length(hits))
	# println(hits)
	for (i, hit) in enumerate(hits)
		filename::String = string(basename, "/", conf, "/", hit)
		fid::HDF5.File = h5open(filename, "r")
		# println(keys(fid)[1],typeof(fid))
		# println(hit)
		for (iTMOS, TMOS) in enumerate(TMOSs)
			group::String = @sprintf("/%s/mesons/%c%.4e_%c%.4e_stoch_stoch", keys(fid)[1], TMOS[1], m, TMOS[2], m)

			raw_data = read(fid, group)
			Gamma_contraction(acc_data, raw_data, T, iTMOS)


		end
		close(fid)
	end

	factor = L^3 * length(hits)
	for (iTMOS, TMOS) in enumerate(TMOSs)
		for t in 1:T
			for i in 1:10
				acc_data[t, i, iTMOS, 1] /= factor
				acc_data[t, i, iTMOS, 2] /= factor
			end
		end
	end
	symm_t!(acc_data)

	println("stoc stoc")
	println("1:  ", acc_data[1, 1, 1, 1], "  +I   ", acc_data[1, 1, 1, 2])
	println("1 T:  ", acc_data[T, 1, 1, 1], "  +I   ", acc_data[T, 1, 1, 2])
	println("1 OS:    ", acc_data[1, 1, 2, 1], "  +I   ", acc_data[1, 1, 2, 2])
	println("1 OS T:  ", acc_data[T, 1, 2, 1], "  +I   ", acc_data[T, 1, 2, 2])
	println("2:  ", acc_data[1:3, 2, 1, 1])
	println("3:  ", acc_data[1:3, 3, 1, 1])


end

struct struct_TM
end
struct struct_OS
end

function eigenvalue_littleD(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_TM)
	return ((a + (b)im) * (c - (d)im))
end
function eigenvalue_littleD(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_OS)
	return ((a + (b)im) * (c + (d)im))
end

struct struct_Esq_her
end
struct struct_Esq_antiher
end

function Esq(a::Float64, b::Float64, c::Float64, d::Float64)
	return ((a + (b)im) * (c - (d)im))
end
function Esq(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_Esq_her)
	return ((a + (b)im) * (c - (d)im))
end
function Esq(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_Esq_antiher)
	return (-(a + (b)im) * (c - (d)im))
end

function convolution_E(ct::Array{ComplexF64, 1}, P::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, Tuple{Int64}}, raw_data::Array{Float64, 5}, ig::Int, count::Int, T::Int)
	for t in 1:T
		ct[t] = raw_data[1, ig, count, 1, t] + (raw_data[2, ig, count, 1, t])im
	end
	cp = P * ct

	# cp[1] = cp[1] * conj(cp[1])
	for p in 1:T
		# cp[p] = cp[p] * conj(cp[p])
		cp[p] = abs2(cp[p])
	end

	return ((P * cp))
end


function exact_exact(conf::String, basename::String, L::Int, T::Int, nvec::Int, TMOSs::Vector{Vector{String}})
	littleD::Array{Float64, 3} = Array{Float64, 3}(undef, 2, nvec, nvec)
	filename = basename * "/" * conf * "/twop_exact_exact_nvec" * string(nvec) * ".h5"
	fid = h5open(filename, "r")
	# println(keys(fid)[1])
	littleD = read(fid, "littleD")

	dimT::Int = div(nvec * nvec - nvec, 2) + nvec
	# {128, 1, 80200, 16, 2}
	raw_data::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, dimT, 1, T)
	raw_data = read(fid, "/sx00sy00sz00st00/mesons/exact_exact")

	# ct::Array{Float64, 2} = zeros(Float64, 2, T)
	ct::Array{ComplexF64, 1} = zeros(ComplexF64, T)
	ct_acc::Array{ComplexF64, 1} = zeros(ComplexF64, T)
	ctave::Array{Float64, 4} = zeros(T, 10, 2, 2)

	count::Int = 1
	P::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, Tuple{Int64}} = plan_fft(ct)

	for (iTMOS, TMOS) in enumerate(["+", "+"]) # enumerate(TMOSs)
		local sign_TMOS#:: Union{struct_TM,struct_OS}
		if (iTMOS == 1)
			sign_TMOS = struct_TM()
		else
			sign_TMOS = struct_OS()
		end
		#"pseudoscalar, scalar, g5g1, g5g2, g5g3, g5g4, g1, g2, g3, g4,s12,s13,s23,s41,s42,s43"
		for ig in 1:10 #1:10
			local is_her#:: Union{struct_Esq_antiher,struct_Esq_her}
			if (ig > 2 && ig < 7)
				is_her = struct_Esq_antiher()
			else
				is_her = struct_Esq_her()
			end
			fill!(ct_acc, 0.0)
			#######################
			## upper part
			#######################
			count = 1
			for i::Int in 1:(nvec-1)
				count += 1
				for j::Int in (i+1):nvec
					ct = convolution_E(ct, P, raw_data, ig, count, T)
					r2::ComplexF64 = eigenvalue_littleD(littleD[1, i, i], littleD[2, i, i], littleD[1, j, j], littleD[2, j, j], sign_TMOS)
					# r2=1
					ct_acc .+= ct / r2
					# ct_acc .+= real(ct / r2)
					count += 1
				end
			end
			for t in 1:T
				ct_acc[t] *= 2
			end
			#######################
			## diag
			#######################
			count = 1
			for i in 1:nvec
				ct = convolution_E(ct, P, raw_data, ig, count, T)
				r2::ComplexF64 = eigenvalue_littleD(littleD[1, i, i], littleD[2, i, i], littleD[1, i, i], littleD[2, i, i], sign_TMOS)
				# r2=1
				ct_acc .+= ct / r2
				count += nvec - i + 1
			end
			factor::Float64 = T * T * L^3 # extra T because we did not put it in the FFT
			for t in 1:T
				ctave[t, ig, iTMOS, 1] = ct_acc[t].re / factor
				ctave[t, ig, iTMOS, 2] = ct_acc[t].im / factor
			end
		end
	end

	symm_t!(ctave)
	# 1,   ,    ,nvec
	#  ,nvec+1, ,2nvec-2
	println("exact_exact  ")
	println("1:  ", ctave[1, 1, 1, 1], "  +I   ", ctave[1, 1, 1, 2])
	println("1 T:  ", ctave[T, 1, 1, 1], "  +I   ", ctave[T, 1, 1, 2])
	println("1 OS:    ", ctave[1, 1, 2, 1], "  +I   ", ctave[1, 1, 2, 2])
	println("1 OS T:  ", ctave[T, 1, 2, 1], "  +I   ", ctave[T, 1, 2, 2])


	println("2:  ", ctave[1, 2, 1, 1], "  +I   ", ctave[1, 2, 1, 2])
	println("2 T:  ", ctave[T, 2, 1, 1], "  +I   ", ctave[T, 2, 1, 2])
	

	# println(raw_data[1,2,1],"  ",raw_data[1,2,2])
	# println(raw_data[2,1,1],"  ",raw_data[2,1,2])
	# println(raw_data[2,2,1],"  ",raw_data[2,2,2])
	# exit(1)
	close(fid)

end

function main()
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/B64/pion_defl/"
	T::Int = 128
	L::Int = 64
	beta::Float64 = 1.778000000000 ##check
	kappa::Float64 = 0.139426500000 ##check
	mass::Float64 = 7.2000e-04
	TMOSs = [["+", "+"], ["+", "-"]]
	nvec::Int = 400

	##find confs
	println(basename)
	confs = readdir(basename)
	conf_new = findall(occursin.(r"^[0-9][0-9][0-9][0-9]_r[0-9]$", confs))
	confs = confs[conf_new]
	println("confs: ", length(confs))

	acc_data::Array{Float64, 4} = zeros(T, 10, 2, 2)
	for (iconf, conf) in enumerate(["2140_r0"])#enumerate(confs)
		println(conf)

		# exact_exact(conf, basename, T, nvec, confs)
		@time exact_exact(conf, basename, L, T, nvec, TMOSs)

		@time stoch_stoch(conf, acc_data, basename, L, T, nvec, mass, TMOSs)
	end

end

main()