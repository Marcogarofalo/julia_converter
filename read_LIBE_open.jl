using HDF5
using Printf
# import Base.write
# using Profile
# using Statistics

# include("./read_hdf5.jl")
# qui i correlatori hanno 5 componenti
# le combinazioni sono 00, 01, 10, 02, 20
# dove 0 e' il propagator normale, 1 e' D^2 e 2 e' D g5 D
# include("./gamma.jl")

# corr::Array{Float64, 7} = zeros(Float64, length(masses), 2, length(TMOSs), length(gammas), length(counterterms), T, 2)

function Gamma_contraction_minimal!(acc_data::Array{Float64, 3}, raw_data::Array{Float64, 6}, T::Int32, g5GI_list::Array{gamma_struct,2})

	# Threads.@threads
	for t in 1:T
		for ig in 1:size(g5GI_list)[1]
			# for ig2 in 1:length(g5GI_list)
				for j in 1:4
					for k in 1:4
						tmp::Complex{Float64} = (raw_data[1, j, g5GI_list[ig,1].col[k], k, g5GI_list[ig,2].col[j], t]
												 +
												 (raw_data[2, j, g5GI_list[ig,1].col[k], k, g5GI_list[ig,2].col[j], t])im) * g5GI_list[ig,2].val[j] * g5GI_list[ig,1].val[k]
						acc_data[ig, t, 1] += tmp.re
						acc_data[ig, t, 2] += tmp.im
					end
				end
			# end
		end
	end
end


function Gamma_contraction(acc_data::Array{Float64, 7}, im::Int, im1::Int, iTMOS::Int, ic::Int, raw_data::Array{Float64, 9}, T::Int32, g5GI_list::Vector{gamma_struct})

	ie1 = (ic - 1) % 3 + 1
	ie2 = div(ic - 1, 3) + 1
	for t in 1:T
		for ig1 in 1:length(g5GI_list)
			for ig2 in 1:length(g5GI_list)
				for j in 1:4
					for k in 1:4
						tmp::Complex{Float64} =
							(raw_data[1, j, g5GI_list[ig1].col[k], k, g5GI_list[ig2].col[j], ie1, ie2, 1, t]
							 +
							 (raw_data[2, j, g5GI_list[ig1].col[k], k, g5GI_list[ig2].col[j], ie1, ie2, 1, t])im) *
							g5GI_list[ig2].val[j] * g5GI_list[ig1].val[k]
						acc_data[im, im1, iTMOS, ig1+(ig2-1)*length(g5GI_list), ic, t, 1] += tmp.re
						acc_data[im, im1, iTMOS, ig1+(ig2-1)*length(g5GI_list), ic, t, 2] += tmp.im
					end
				end
			end
		end
	end
end

function copy_slice(Ng::Int32, im::Int32, im1::Int32, iTMOS::Int32, ic::Int32)

	# for i in length(g5G)
	# 	for t in 1:T
	# 		slice[i,t,2]= acc_data[im, im1, iTMOS, i, ic, 1, t]
end


function read_LIBE_open(conf::String, acc_data::Array{Float64, 7}, basename::String, L::Int32, T::Int32, masses::Vector{Float64},
	TMOSs::Vector{Vector{String}}, info_counterterms::Vector{Int32}, g5GI_list::Array{gamma_struct,2})
	# {96, 1, 3, 3, 4, 4, 4, 4, 2}
	# mu1_mu2_{       e1, e2} --- > so in julia they become as {e2,e1}
	raw_data::Array{Float64, 9} = Array{Float64, 9}(undef, 2, 4, 4, 4, 4, 3, 3, 1, T)
	slice_raw::Array{Float64, 6} = Array{Float64, 6}(undef, 2, 4, 4, 4, 4, T)
	slice::Array{Float64, 3} = Array{Float64, 3}(undef, size(g5GI_list)[1], T, 2)

	#sib 
	# { 96, 1, 5, 4, 4, 4, 4, 2}
	raw_sib::Array{Float64, 8} = Array{Float64, 8}(undef, 2, 4, 4, 4, 4, 5, 1, T)


	hits::Vector{String} = readdir(string(basename, "/", conf))
	pattern::String = string("^twop_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new]
	println("Nhits:   ", length(hits))
	m1list::Array{Float64, 1} = [masses[1], masses[1]]

	for (i, hit) in enumerate(hits)
		filename::String = string(basename, "/", conf, "/", hit)
		fid::HDF5.File = h5open(filename, "r")
		for (im, m) in enumerate(masses)
			m1list[2] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)

					group::String = @sprintf("/%s/mesons/%c%.4e_%c%.4e_qed_open", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					raw_data = read(fid, group)

					for ic in 1:info_counterterms[1]
						# Gamma_contraction(acc_data, im, im1, iTMOS, ic, raw_data, T, g5GI_list)
						ie1::Int = (ic - 1) % 3 + 1
						ie2::Int = div(ic - 1, 3) + 1
						slice_raw .= raw_data[:, :, :, :, :, ie2, ie1, 1, :]
						fill!(slice, 0.0)#slice= acc_data[im, im1, iTMOS, :, ic, :, :]
						Gamma_contraction_minimal!(slice, slice_raw, T, g5GI_list)
						acc_data[im, im1, iTMOS, :, ic, :, :] += slice
					end
					
					group = @sprintf("/%s/mesons/%c%.4e_%c%.4e_sib_open", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					raw_sib = read(fid, group)

					for ic in (info_counterterms[1]+1):(info_counterterms[1]+info_counterterms[2])
						# Gamma_contraction(acc_data, im, im1, iTMOS, ic, raw_data, T, g5GI_list)
						is::Int = ic - info_counterterms[1]
						slice_raw .= raw_sib[:, :, :, :, :, is, 1, :]
						slice .= acc_data[im, im1, iTMOS, :, ic, :, :]
						Gamma_contraction_minimal!(slice, slice_raw, T, g5GI_list)
						acc_data[im, im1, iTMOS, :, ic, :, :] .= slice
					end


				end
			end


		end
		close(fid)
	end

	factor::Float64 = L^3 * length(hits)
	acc_data /= factor
	# for (iTMOS, TMOS) in enumerate(TMOSs)
	# 	for t in 1:T
	# 		for i in 1:16
	# 			for j in 1:16
	# 				for ins in 1:5
	# 					acc_data[t, i, j, ins, iTMOS, 1] /= factor
	# 					acc_data[t, i, j, ins, iTMOS, 2] /= factor
	# 				end
	# 			end
	# 		end
	# 	end
	# end
	# symm_t!(acc_data)

	# println("stoc stoc")
	# for i in 1:16
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 0, acc_data[1, i, 1, 1], acc_data[1, i, 1, 2], acc_data[1, i, 2, 1], acc_data[1, i, 2, 2])
	# 	@printf("id:%-3d  t:%-3d  TM   %-20.12g  +I %-20.12g   OS   %-20.12g  +I %-20.12g\n", (i - 1), 1, acc_data[2, i, 1, 1], acc_data[2, i, 1, 2], acc_data[2, i, 2, 1], acc_data[2, i, 2, 2])
	# end
	GC.gc() 

end

function read_LIBE!(conf::String, acc_data::Array{Float64, 8}, iconf::Int, basename::String, L::Int32, T::Int32, masses::Vector{Float64},
	TMOSs::Vector{Vector{String}}, info_counterterms::Vector{Int32}, g5GI_list::Vector{gamma_struct})
	# {96, 1, 3, 3, 16, 2}
	# mu1_mu2_{       e1, e2} --- > so in julia they become as {e2,e1}
	# raw_data::Array{Float64, 6} = Array{Float64, 6}(undef, 2, length(g5GI_list), 3, 3, 1, T)
	# slice_raw::Array{Float64, 6} = Array{Float64, 6}(undef, 2, 16, T)
	# slice::Array{Float64, 3} = Array{Float64, 3}(undef, length(g5GI_list), T, 2)
	raw_data::Array{Float64, 6} = Array{Float64, 6}(undef, 2, length(g5GI_list), 3, 3, 1, T)

	#sib 
	# {96, 1, 5, 16, 2}
	raw_sib::Array{Float64, 5} = Array{Float64, 5}(undef, 2, 16, 5, 1, T)


	hits::Vector{String} = readdir(string(basename, "/", conf))
	pattern::String = string("^twop_id[0-9]*_st[0-9]*\\.h5\$")
	local conf_new = findall(occursin.(Regex(pattern), hits))
	hits = hits[conf_new]
	println(conf,"  Nhits:   ", length(hits))
	m1list::Array{Float64, 1} = [masses[1], masses[1]]

	for (i, hit) in enumerate(hits)
		filename::String = string(basename, "/", conf, "/", hit)
		fid::HDF5.File = h5open(filename, "r")
		for (im, m) in enumerate(masses)
			m1list[2] = m
			for (im1, m1) in enumerate(m1list)
				for (iTMOS, TMOS) in enumerate(TMOSs)

					group::String = @sprintf("/%s/mesons/%c%.4e_%c%.4e_qed", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					raw_data = read(fid, group)
					
					for ic in 1:info_counterterms[1]
						# Gamma_contraction(acc_data, im, im1, iTMOS, ic, raw_data, T, g5GI_list)
						ie1::Int = (ic - 1) % 3 + 1
						ie2::Int = div(ic - 1, 3) + 1
						for t in 1:T
							for ig in 1:length(g5GI_list)
								for reim in 1:2
									acc_data[iconf, im, im1, iTMOS, ig, ic, t, reim] += raw_data[reim, ig, ie2, ie1, 1, t]
								end
							end
						end
					end

					group = @sprintf("/%s/mesons/%c%.4e_%c%.4e_sib", keys(fid)[1], TMOS[1], m, TMOS[2], m1)
					raw_sib = read(fid, group)

					for ic in (info_counterterms[1]+1):(info_counterterms[1]+info_counterterms[2])
						# Gamma_contraction(acc_data, im, im1, iTMOS, ic, raw_data, T, g5GI_list)
						is::Int = ic - info_counterterms[1]
						# slice_raw .= raw_sib[:, :, :, :, :, is, 1, :]
						# slice .= acc_data[im, im1, iTMOS, :, ic, :, :]
						# Gamma_contraction_minimal!(slice, slice_raw, T, g5GI_list)
						# acc_data[im, im1, iTMOS, :, ic, :, :] .= slice
						for t in 1:T
							for ig in 1:length(g5GI_list)
								for reim in 1:2
									acc_data[iconf, im, im1, iTMOS, ig, ic, t, reim] += raw_sib[reim, ig, is, 1, t]
								end
							end
						end
					end
					
				end
			end


		end
		close(fid)
	end

	factor::Float64 = L^3 * length(hits)
	# acc_data[iconf, : ] /= factor
	for (im, m) in enumerate(masses)
		m1list[2] = m
		for (im1, m1) in enumerate(m1list)
			for (iTMOS, TMOS) in enumerate(TMOSs)
				for ic in 1:(info_counterterms[1]+info_counterterms[2])
					for t in 1:T
						for ig in 1:length(g5GI_list)
							for reim in 1:2
								acc_data[iconf, im, im1, iTMOS, ig, ic, t, reim] /= factor
							end
						end
					end
				end
				if (TMOS[1]=="+" && TMOS[2]=="+" && m==m1 && im==length(masses))
					ic = 4+1
					println("data 4 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					ic = 9+1
					println("data 9 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					ic = 10 +1
					println("data 10 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					ic = 11 +1
					println("data 11 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					ic = 12 +1
					println("data 12 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					ic = 13 +1
					println("data 13 ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 1], " ", acc_data[iconf, im, im1, iTMOS, 1, ic, 1, 2]  )
					
					println()
				end
			end
		end
	end
	GC.gc() 
end
