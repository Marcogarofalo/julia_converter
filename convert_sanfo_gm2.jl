#print(map(x->string(x, x), ARGS))
using Printf

if length(ARGS)!=2
    println("usage: julia convert_sanfo_gmw.jl  binfile  gm2file")
    exit(1)
end
if ARGS[1]==ARGS[2]
    println("input and output must be different")
    exit(1)
end

println("reading: ",ARGS[1])
fs=filesize(ARGS[1])
println("size: ",fs)
io = open(ARGS[1], "r");

Nconfs= read(io, Int64)
Nhits = read(io, Int64)
T     = read(io, Int64)
println("Nconfs: ",Nconfs)
println("Nhits: ",Nhits)
println("T: ",T)

if (fs-sizeof(Int64)*3-(Nconfs*sizeof(Float64)*(T/2+1)) !=0)
    println("size and Nconfs do not match")
    exit(1)
end

chunk_size=Int(T/2+1)
#data = zeros(Float64, (Nconfs,chunk_size ))
# data=Array{Float32}(undef, chunk_size, Nconfs)
data=Array{Float64}(undef, Nconfs, chunk_size)

for i in 1:Nconfs
    for t in 1:chunk_size    
    data[i,t]=read(io, Float64)
    end
end
# wrong order to read it in one go
# read!(io, data)

close(io)

out = open(ARGS[2], "w");
for i in 1:Nconfs
    println(out,"")
    @printf(out,"# %04d_r0\n",i)
    println(out,"")
    for t in 1:chunk_size
        println(out,data[i,t])
        # println(out,data[i,t])
    end
end
   

