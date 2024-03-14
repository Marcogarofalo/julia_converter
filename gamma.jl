
const id::Matrix{Complex{Float64}} = [
	1.0 0.0 0.0 0.0;
	0.0 1.0 0.0 0.0;
	0.0 0.0 1.0 0.0;
	0.0 0.0 0.0 1.0]

const gx::Matrix{Complex{Float64}} = [
	 0  0 0 im;
	 0  0 im 0;
	 0 -im 0 0;
	-im  0 0 0]

const gy::Matrix{Complex{Float64}} = [
	0  0  0  1;
	0  0 -1  0;
	0 -1  0  0;
	1  0  0  0]


const gz::Matrix{Complex{Float64}} = [
	 0  0  im  0;
	 0  0  0  -im;
	-im 0  0   0;
	 0  im 0   0]

const gt::Matrix{Complex{Float64}} = [
	1  0  0  0;
	0  1  0  0;
	0  0 -1  0;
	0  0  0 -1]

const g5::Matrix{Complex{Float64}} = [
	 0  0 -1  0;
	 0  0  0 -1;
	-1  0  0  0;
	 0 -1  0  0]

# const gamma = [id, g5, gx, gy, gz, gt, g5 * gx, g5 * gy, g5 * gz, g5 * gt]
# this will go in p5,id, p5Vk, VK

struct gamma_struct
	col::Vector{Int32}
	val::Vector{ComplexF64}
end

const G5::gamma_struct=gamma_struct(
	 [3, 4, 1, 2],
	 [-1, -1, -1, -1]
)

const Gt::gamma_struct=gamma_struct(
	 [1, 2, 3, 4],
	 [1, 1, -1, -1]
	)

const Gz::gamma_struct=gamma_struct(
	 [3, 4, 1, 2],
	 [im, -im, -im, im]
	)

const Gy::gamma_struct=gamma_struct(
	 [4, 3, 2, 1],
	 [1, -1, -1, 1]
	)
const Gx::gamma_struct=gamma_struct(
	 [4, 3, 2, 1],
	 [im, im, -im, -im]
	)
const Id::gamma_struct=gamma_struct(
	 [1,2,3,4],
	 [1, 1,1,1]
)

function get_matrix_form(a::gamma_struct)
	Ma=zeros(ComplexF64, 4, 4)
	for i in 1:4
		Ma[i,a.col[i]]=a.val[i]
	end
	return(Ma)
end

function get_struct_form(Mab::Matrix{ComplexF64})
	col=zeros(Int64,4)
	val=zeros(ComplexF64,4)
	for i in 1:4
		for j in 1:4
			if (abs(Mab[i,j])>1e-7)
				col[i]=j
				val[i]=Mab[i,j]
            end
		end
	end
	return(gamma_struct(col,val))
end

function Base.:*(a::gamma_struct, b::gamma_struct) 
	Ma=get_matrix_form(a)
	Mb=get_matrix_form(b)
	
	Mab=Ma*Mb
	
	return(get_struct_form(Mab))
end

function Base.:*(a::Number, b::gamma_struct) 
	Mb=get_matrix_form(b)
	
	Mab=Mb*a
	
	return(get_struct_form(Mab))
end
function Base.:*( b::gamma_struct,a::Number) 
	Mb=get_matrix_form(b)
	
	Mab=Mb*a
	
	return(get_struct_form(Mab))
end
function Base.:+(a::gamma_struct, b::gamma_struct) 
	Ma=get_matrix_form(a)
	Mb=get_matrix_form(b)
	
	Mab=Ma+Mb
	
	return(get_struct_form(Mab))
end

function Base.:-(a::gamma_struct, b::gamma_struct) 
	Ma=get_matrix_form(a)
	Mb=get_matrix_form(b)
	
	Mab=Ma-Mb
	
	return(get_struct_form(Mab))
end

const sxy::gamma_struct= (-im/2) * (Gx*Gy-Gy*Gx)
const sxz::gamma_struct= (-im/2) * (Gx*Gz-Gz*Gx)
const syz::gamma_struct= (-im/2) * (Gy*Gz-Gz*Gy)
const stx::gamma_struct= (-im/2) * (Gt*Gx-Gx*Gt)
const sty::gamma_struct= (-im/2) * (Gt*Gy-Gy*Gt)
const stz::gamma_struct= (-im/2) * (Gt*Gz-Gz*Gt)

# this is g5*Gamma of the LIBE
# sthe stoch-stoch part does not have the g5 inside 
const Gamma::Vector{gamma_struct} = [Id, Gx, Gy, Gz, Gt, G5,
			G5*Gx, G5*Gy, G5*Gz, G5*Gt
			 , sxy, sxz, syz, stx, im*sty, stz]


# const Gamma = [Id, G5, Gx, Gy, Gz, Gt, G5 * Gx, G5 * Gy, G5 * Gz, G5 * Gt,
#   Id,Id,Id,Id,Id,Id]
