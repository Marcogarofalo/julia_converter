#print(map(x->string(x, x), ARGS))
using Printf

if length(ARGS) != 2
	println("usage: julia create_fake_gm2.jl T outfile")
	exit(1)
end

# println("reading: ",ARGS[1])
# fs=filesize(ARGS[1])
# println("size: ",fs)
# io = open(ARGS[1], "r");

Nconfs = 50
Nhits  = 1
T      = parse(Int, ARGS[1])
println("Nconfs: ", Nconfs)
println("Nhits: ", Nhits)
println("T: ", T)
fs = sizeof(Int64) * 3 + (Nconfs * sizeof(Float64) * (T / 2 + 1))

if (fs - sizeof(Int64) * 3 - (Nconfs * sizeof(Float64) * (T / 2 + 1)) != 0)
	println("size ",fs," and Nconfs do not match   sizeof(Int64) * 3 - (Nconfs * sizeof(Float64) * (T / 2 + 1))=",sizeof(Int64) * 3 - (Nconfs * sizeof(Float64) * (T / 2 + 1)))
	exit(1)
end

chunk_size = Int(T / 2 + 1)
#data = zeros(Float64, (Nconfs,chunk_size ))
# data=Array{Float32}(undef, chunk_size, Nconfs)
data = Array{Float64}(undef, Nconfs, chunk_size)

for i in 1:Nconfs
	for t in 1:chunk_size
		data[i, t] = 0 #read(io, Float64)
	end
end
# wrong order to read it in one go
# read!(io, data)

# close(io)

out = open(ARGS[2], "w");
for i in 1:Nconfs
	println(out, "")
	@printf(out, "# %04d_r0\n", i)
	println(out, "")
	for t in 1:chunk_size
		println(out, data[i, t])
		# println(out,data[i,t])
	end
end


