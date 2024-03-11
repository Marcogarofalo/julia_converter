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
			for i in 1:16
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
		# for i in 1:16
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
		for i in 1:16
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

	factor::Float64 = L^3 * length(hits)
	for (iTMOS, TMOS) in enumerate(TMOSs)
		for t in 1:T
			for i in 1:16
				acc_data[t, i, iTMOS, 1] /= factor
				acc_data[t, i, iTMOS, 2] /= factor
			end
		end
	end
	symm_t!(acc_data)

	# println("stoc stoc")
	# for i in 1:16
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 0, acc_data[1, i, 1, 1], acc_data[1, i, 1, 2], acc_data[1, i, 2, 1], acc_data[1, i, 2, 2])
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 1, acc_data[2, i, 1, 1], acc_data[2, i, 1, 2], acc_data[2, i, 2, 1], acc_data[2, i, 2, 2])
	# end
end

struct struct_TM
end
struct struct_OS
end

function inv_eigenvalue_littleD_aver(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_TM)
	r2::ComplexF64 = ((a + (b)im) * (c - (d)im))
	r2 = 1 / r2
	# r2 += conj(r2)
	# r2 /= 2.0
	return (r2.re)
end
function inv_eigenvalue_littleD_aver(a::Float64, b::Float64, c::Float64, d::Float64, type::struct_OS)
	r2::ComplexF64 = ((a + (b)im) * (c + (d)im))
	r2 = 1 / r2
	# r2 += conj(r2)
	# r2 /= 2.0
	return (r2.re)
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
	# ct = P * cp
	return ((P * cp))
end


function exact_exact(conf::String, ctave::Array{Float64, 4}, basename::String, L::Int, T::Int, nvec::Int, TMOSs::Vector{Vector{String}})
	littleD::Array{Float64, 3} = Array{Float64, 3}(undef, 2, nvec, nvec)
	filename = basename * "/" * conf * "/twop_exact_exact_nvec" * string(nvec) * ".h5"
	fid = h5open(filename, "r")
	# println(keys(fid)[1])
	littleD = read(fid, "littleD")

	lambda_inv::Array{ComplexF64, 1} = Array{Float64, 1}(undef, nvec)
	for i::Int in 1:nvec
		lambda_inv[i] = 1.0 / (littleD[1, i, i] + littleD[2, i, i]im)
	end

	dimT::Int = div(nvec * nvec - nvec, 2) + nvec
	# {128, 1, 80200, 16, 2}
	raw_data::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, dimT, 1, T)
	raw_data = read(fid, "/sx00sy00sz00st00/mesons/exact_exact")

	# ct::Array{Float64, 2} = zeros(Float64, 2, T)
	ct::Array{ComplexF64, 1} = zeros(ComplexF64, T)
	cp::Array{ComplexF64, 1} = zeros(ComplexF64, T)
	ct_acc::Array{ComplexF64, 2} = zeros(ComplexF64, T, 2)
	# ctave::Array{Float64, 4} = zeros(T, 16, 2, 2)

	count::Int = 1
	P::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, Tuple{Int64}} = plan_fft(ct)

	# for (iTMOS, TMOS) in enumerate(TMOSs)
	# 	local sign_TMOS#:: Union{struct_TM,struct_OS}
	# 	if (iTMOS == 1)
	# 		sign_TMOS = struct_TM()
	# 	else
	# 		sign_TMOS = struct_OS()
	# 	end
	#"pseudoscalar, scalar, g5g1, g5g2, g5g3, g5g4, g1, g2, g3, g4,s12,s13,s23,s41,s42,s43"
	for ig in 1:16 #1:16
		fill!(ct_acc, 0.0)
		#######################
		## upper part
		#######################
		count = 1
		for i::Int in 1:(nvec-1)
			count += 1
			for j::Int in (i+1):nvec
				ct = convolution_E(ct, P, raw_data, ig, count, T)
				# r2::Float64 = inv_eigenvalue_littleD_aver(littleD[1, i, i], littleD[2, i, i], littleD[1, j, j], littleD[2, j, j], sign_TMOS)
				r2::ComplexF64 = (lambda_inv[i] * conj(lambda_inv[j]))
				for t in 1:T
					# if (abs(ct[t].im)>0)
					# 	println("  ERROR    ",t," ",ct[t]*r2)
					# end
					#there is no need to sum the other r combination, it will be exactly the same
					ct_acc[t, 1] += ct[t] * (r2 + conj(r2))

				end
				r2 = (lambda_inv[i] * lambda_inv[j])
				for t in 1:T
					# the other combination i<->j is the same, so x2 but then we do the average with +-r  
					ct_acc[t, 2] += ct[t] * (r2 + conj(r2))
					# r=-1-1
					# ct_acc[t,2] += 2*ct[t] * conj(r2)
					# ct_acc[t,2] /= 2

				end
				# ct_acc .+= ct * r2
				count += 1
			end
		end
		# for t in 1:T
		# 	ct_acc[t,1] *= 2
		# 	ct_acc[t,2] *= 2
		# end
		#######################
		## diag
		#######################
		count = 1
		for i in 1:nvec
			ct = convolution_E(ct, P, raw_data, ig, count, T)

			# r2::Float64 = inv_eigenvalue_littleD_aver(littleD[1, i, i], littleD[2, i, i], littleD[1, i, i], littleD[2, i, i], sign_TMOS)
			#TM
			r2::ComplexF64 = (lambda_inv[i] * conj(lambda_inv[i]))
			for t in 1:T
				# if (abs(ct[t].im)>0)
				# 	println("  ERROR diag   ",t," ",ct[t]*r2)
				# end
				ct_acc[t, 1] += ct[t] * r2
			end
			#OS
			r2 = (lambda_inv[i] * lambda_inv[i])
			for t in 1:T
				ct_acc[t, 2] += ct[t] * (r2 + conj(r2)) / 2.0
			end
			count += nvec - i + 1
		end
		factor::Float64 = T * T * L^3 # extra T because we did not put it in the FFT
		if (ig >= 7 && ig <= 10)
			factor *= -1
		end
		for t in 1:T
			#TM
			ctave[t, ig, 1, 1] = ct_acc[t, 1].re / factor
			ctave[t, ig, 1, 2] = ct_acc[t, 1].im / factor
			#OS
			ctave[t, ig, 2, 1] = ct_acc[t, 2].re / factor
			ctave[t, ig, 2, 2] = ct_acc[t, 2].im / factor
		end

		# end
	end

	symm_t!(ctave)
	# 1,   ,    ,nvec
	#  ,nvec+1, ,2nvec-2
	# println("exact_exact  ")
	# for i in 1:16
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 0, ctave[1, i, 1, 1], ctave[1, i, 1, 2], ctave[1, i, 2, 1], ctave[1, i, 2, 2])
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 1, ctave[2, i, 1, 1], ctave[2, i, 1, 2], ctave[2, i, 2, 1], ctave[2, i, 2, 2])
	# end
	# println(raw_data[1,2,1],"  ",raw_data[1,2,2])
	# println(raw_data[2,1,1],"  ",raw_data[2,1,2])
	# println(raw_data[2,2,1],"  ",raw_data[2,2,2])
	# exit(1)
	close(fid)

end

function accumulate_tvg!(ctave::Array{Float64, 4}, littleD::Array{Float64, 3},
	data_p::Array{Float64, 5}, data_m::Array{Float64, 5})

	for t in 1:size(ctave)[1]
		for iv in 1:size(littleD)[2]
			rp::ComplexF64 = littleD[1, iv, iv] + littleD[2, iv, iv]im
			rm::ComplexF64 = littleD[1, iv, iv] - littleD[2, iv, iv]im
			for ig in 1:16
				dp::Float64 = data_p[1, ig, iv, 1, t] #+ data_p[2, ig, iv, 1, t]im
				dm::Float64 = data_m[1, ig, iv, 1, t] #+ data_m[2, ig, iv, 1, t]im
				# ee::ComplexF64 = data_e[1,ig,iv,1,t]
				# println(dp,"  ",dm,"         ",rp,"   ",rm)
				rOS::ComplexF64 = dp / rm + dm / rp
				# rOS::ComplexF64 = dp / rm + conj(dm) / rp + conj(dp) / rm + dm / rp;
				# rTM::ComplexF64 = dm / rm + dp / rp
				# rTM::ComplexF64 = dm / rm + dp / rp + conj(dm) / rm + conj(dp) / rp ;
				ctave[t, ig, 1, 1] += rOS.re#rTM.re
				ctave[t, ig, 1, 2] += -rOS.im#rTM.im
				ctave[t, ig, 2, 1] += rOS.re
				ctave[t, ig, 2, 2] += rOS.im
			end
		end
	end

end

function mult_eigen!(ctave::Array{Float64, 4}, littleD::Array{Float64, 3},
	data_p_ave::Array{ComplexF64, 3}, data_m_ave::Array{ComplexF64, 3})

	for t in 1:size(ctave)[1]
		for iv in 1:size(littleD)[2]
			rp::ComplexF64 = littleD[1, iv, iv] + littleD[2, iv, iv]im
			rm::ComplexF64 = littleD[1, iv, iv] - littleD[2, iv, iv]im
			for ig in 1:16
				dp::ComplexF64 = data_p_ave[ig, iv, t]
				dm::ComplexF64 = data_m_ave[ig, iv, t]

				rOS::Float64 = real(dp / rm) + real(dm / rp)
				rTM::Float64 = real(dm / rm) + real(dp / rp)

				ctave[t, ig, 1, 1] += rTM
				# ctave[t, ig, 1, 2] += rTM.im
				ctave[t, ig, 2, 1] += rOS
				# ctave[t, ig, 2, 2] += rOS.im
			end
		end
	end

end

function stoch_exact(conf::String, ctave::Array{Float64, 4}, basename::String, L::Int, T::Int, nvec::Int, m::Float64, TMOSs::Vector{Vector{String}})
	littleD::Array{Float64, 3} = Array{Float64, 3}(undef, 2, nvec, nvec)
	filename_exact::String = basename * "/" * conf * "/twop_exact_exact_nvec" * string(nvec) * ".h5"
	fid = h5open(filename_exact, "r")
	# println(keys(fid)[1])
	littleD .= read(fid, "littleD")
	close(fid)
	#{128, 1, 400, 16, 2}
	data_p::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, nvec, 1, T)
	data_m::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, nvec, 1, T)

	data_p_ave::Array{ComplexF64, 3} = zeros(ComplexF64, 16, nvec, T)
	data_m_ave::Array{ComplexF64, 3} = zeros(ComplexF64, 16, nvec, T)

	# ctave::Array{Float64, 4} = zeros(T, 16, 2, 2)

	hits = readdir(string(basename, "/", conf))
	pattern = string("^twop_nev", nvec, "_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new]
	println("Nhits:   ", length(hits))

	# for (i, hit) in enumerate(hits)
	# 	filename::String = string(basename, "/", conf, "/", hit)
	# 	fid_s::HDF5.File = h5open(filename, "r")

	# 	group_p::String = @sprintf("/%s/mesons/+%.4e_stoch_exact_G_G", keys(fid_s)[1], m)
	# 	group_m::String = @sprintf("/%s/mesons/-%.4e_stoch_exact_G_G", keys(fid_s)[1], m)

	# 	data_p .= read(fid_s, group_p)	
	# 	data_m .= read(fid_s, group_m)
	# 	close(fid_s)

	# 	accumulate_tvg!(ctave, littleD, data_p, data_m )

	# end
	for (i, hit) in enumerate(hits)
		filename::String = string(basename, "/", conf, "/", hit)
		fid_s::HDF5.File = h5open(filename, "r")

		group_p::String = @sprintf("/%s/mesons/+%.4e_stoch_exact_G_G", keys(fid_s)[1], m)
		group_m::String = @sprintf("/%s/mesons/-%.4e_stoch_exact_G_G", keys(fid_s)[1], m)

		data_p .= read(fid_s, group_p)
		data_m .= read(fid_s, group_m)
		close(fid_s)
		for t in 1:size(ctave)[1]
			for iv in 1:size(littleD)[2]
				for ig in 1:16

					data_p_ave[ig, iv, t] += data_p[1, ig, iv, 1, t] + data_p[2, ig, iv, 1, t]im
					data_m_ave[ig, iv, t] += data_m[1, ig, iv, 1, t] + data_m[2, ig, iv, 1, t]im
				end
			end
		end
	end

	mult_eigen!(ctave, littleD, data_p_ave, data_m_ave)

	factor::Float64 = L^3 * length(hits)
	for (iTMOS, TMOS) in enumerate(TMOSs)
		for t in 1:size(ctave)[1]
			for i in 1:16
				ctave[t, i, iTMOS, 1] /= factor
				ctave[t, i, iTMOS, 2] /= factor
			end
		end
	end
	symm_t!(ctave)

	# println("stoc_exact  ")
	# for i in 1:16
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 0, ctave[1, i, 1, 1], ctave[1, i, 1, 2], ctave[1, i, 2, 1], ctave[1, i, 2, 2])
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 1, ctave[2, i, 1, 1], ctave[2, i, 1, 2], ctave[2, i, 2, 1], ctave[2, i, 2, 2])
	# end
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

	data_ee::Array{Float64, 4} = zeros(T, 16, 2, 2)
	data_se::Array{Float64, 4} = zeros(T, 16, 2, 2)
	data_ss::Array{Float64, 4} = zeros(T, 16, 2, 2)
	data::Array{Float64, 5} = zeros(length(confs), T, 16, 2, 2)
	for (iconf, conf) in enumerate(confs) #= enumerate(["2140_r0"]) =#
		println(conf)
		fill!(data_ee, 0.0)
		fill!(data_se, 0.0)
		fill!(data_ss, 0.0)

		@time exact_exact(conf, data_ee, basename, L, T, nvec, TMOSs)
		@time stoch_exact(conf, data_se, basename, L, T, nvec, mass, TMOSs)
		@time stoch_stoch(conf, data_ss, basename, L, T, nvec, mass, TMOSs)
		for t in 1:T
			for ig in 1:16
				for iTMOS in 1:2
					for reim in 1:2
						data[iconf, t, ig, iTMOS, reim] = data_ee[t, ig, iTMOS, reim] + data_se[t, ig, iTMOS, reim] + data_ss[t, ig, iTMOS, reim]
					end
				end
			end
		end
	end

end

main()
