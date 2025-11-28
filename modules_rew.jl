module modules_rew

using Printf
using Statistics: Statistics

export char_before_match, compute_wU, monomial_type, TM_monomial, OS_monomial


function char_before_match(s::String, pattern, shift::Integer = Int(0), len_rep_name::Integer = Int(0))
    # println("char_before_match: ", s, " ", pattern, " ", shift, " ", len_rep_name)
	pos = findfirst(pattern, s)
	pos1 = pos[1]
	mys = Vector{String}(undef, 2)  # Define a vector of strings with two elements

	r::Int32 = 0 + shift
	if pos == nothing || pos[1] == 1
		error("Pattern not found or found at the beginning of the string")
	else
		rep = s[(pos1-1):(pos1-1+len_rep_name)]
		# mys[1] = replace(s[(pos1+length(pattern)):length(s)], Regex("onlinemeas.s....") => "")
        mys[1] = s[(length(s)-4):length(s)]
		mys[1] = @sprintf("%04d", parse(Int32, mys[1]))
		if rep == "a"  || rep == "a_cpu"
			mys[2] = mys[1] * "_r" * string(r)
		elseif rep == "b" || rep == "b_cpu"
			mys[2] = mys[1] * "_r" * string(r + 1)
		elseif rep == "c" || rep == "a_gpu"
			mys[2] = mys[1] * "_r" * string(r + 2)
		elseif rep == "d" || rep == "b_gpu"
			mys[2] = mys[1] * "_r" * string(r + 3)
		elseif rep == "e" 
			mys[2] = mys[1] * "_r" * string(r + 4)
		elseif rep == "f" 
			mys[2] = mys[1] * "_r" * string(r + 5)	
		else
			error("Pattern not found or found at the beginning of the string")
		end

		return mys
	end
end

struct TM_monomial end
struct OS_monomial end
monomial_type = Union{TM_monomial, OS_monomial}

function compute_wU(ws::Vector{Float64}, monomial::TM_monomial)
	# println("TM")
	m = Statistics.mean(ws)
	wU = m + log(sum(exp.(ws .- m)))
	return wU
end
function compute_wU(ws::Vector{Float64}, monomial::OS_monomial)
	# println("OS")
	m = Statistics.mean(ws)
	# wU = m / 2.0 + log(sum(exp.((ws .- m) ./ 2.0))) # wrong formula
	wU = (m + log(sum(exp.(ws .- m)))) / 2.0
	return wU
end

end # module modules_rew
