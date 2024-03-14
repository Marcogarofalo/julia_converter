include("corr_builder.jl")

function main()
	basename::String = "/leonardo_scratch/large/userexternal/sbacchio/C80/pion_defl"

	T::Int = 160
	L::Int = 80
	beta::Float64 = 1.778000000000 ##check
	kappa::Float64 = 0.139426500000 ##check
	mass::Float64 = 6.0000e-04
	TMOSs = [["+", "+"], ["+", "-"]]
	nvec::Int = 430

	info::builder_info =builder_info(basename,T,L,beta,kappa,mass,TMOSs,nvec)

	builder(info,"twop","0005_r1")

end

main()
