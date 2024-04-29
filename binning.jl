

function bin_intoN(data::Vector{Float64}, Nb::Int)

	Nconf_in::Float64 = length(data)

	clustSize::Float64 = Nconf_in / Nb
	# // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
	to_write::Vector{Float64} = zeros(Float64, Nb)
	for iClust in 1:Nb

		# /// Initial time of the bin
		binBegin::Float64 = (iClust - 1) * clustSize
		# /// Final time of the bin
		binEnd::Float64 = binBegin + clustSize
		binPos::Float64 = binBegin
		while (binEnd - binPos > 1e-10)
			# /// Index of the configuration related to the time
			iConf::Int = floor(binPos + 1e-10)

			# ///Rectangle left point
			beg::Float64 = binPos

			# /// Rectangle right point
			minend::Float64 = min(binEnd, iConf + 1.0)

			# /// Rectangle horizontal size
			weight::Float64 = minend - beg

			# // Perform the operation passing the info
			to_write[iClust] += weight * data[iConf+1]

			binPos = minend
		end
		to_write[iClust] /= clustSize

	end
	return to_write
end
