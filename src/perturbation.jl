"""
Pertubation switchs binary variables: a switch means that a variable at value 0
will be switched at value 1 and a variable at value 1 will be switched at
value 0.
"""
function noperturb(opts...) end


"""
    randomperturb!(x_tilde::Vector{Float64}, _x_overline::Vector{Float64},
        indices::Vector{Int})

The random perturbation randomly switch variables. Each one has 50% chance to
be switched.
"""
function randomperturb!(
        x_tilde::Vector{Float64},
        _x_overline::Vector{Float64},
        indices::Vector{Int},
        opts...
    )::Vector{Float64}
    for i in indices
        if bitrand()[1]
            x_tilde[i] = abs(x_tilde[i]-1)
        end
    end
    return x_tilde
end


"""
    perturb!(x_tilde::Vector{Float64}, x_overline::Vector{Float64},
        indices::Vector{Int}, TTmin::Int, TTmax::Int)

The perturbation switch the TT most fractional variables. TT is randomly picked
between TTmin to TTmax.
"""
function perturb!(
        x_tilde::Vector{Float64},
        x_overline::Vector{Float64},
        indices::Vector{Int},
        TTmin::Int,
        TTmax::Int,
        opts...
    )::Vector{Float64}
    TT = rand(TTmin:TTmax)

    x_tri = Vector{Tuple{Int, Float64}}()
    for i in indices
        push!(x_tri, (i, abs(x_overline[i] - x_tilde[i])))
    end
    sort!(x_tri, by = x -> x[2], rev=true)

    for i in x_tri[1:TT]
        x_tilde[i[1]] = abs(x_tilde[i[1]]-1)
    end
    return x_tilde
end


"""
    frequencyperturb!(x_tilde::Vector{Float64}, x_overline::Vector{Float64},
        indices::Vector{Int}, TTmin::Int, TTmax::Int, frequency_round::Vector{Int},
        nb_iteration::Int)

The frequency perturbation switch TT variables. The variables are chosen by
frequency and TT is randomly picked between TTmin to TTmax..
"""
function frequencyperturb!(
        x_tilde::Vector{Float64},
        x_overline::Vector{Float64},
        indices::Vector{Int},
        TTmin::Int,
        TTmax::Int,
        frequency_round::Vector{Int},
        nb_iteration::Int
    )::Vector{Float64}
    TT = rand(TTmin:TTmax)

    x_tri = Vector{Tuple{Int, Float64}}()
    for i in indices
        if x_tilde[i] < 0.5
            push!(x_tri, (i, nb_iteration - freq_round[i]))
        else
            push!(x_tri, (i, frequency_round[i]))
        end

    end
    sort!(x_tri, by = x -> x[2], rev=true)

    for i in x_tri[1:TT]
        x_tilde[i[1]] = abs(x_tilde[i[1]]-1)
    end
    return x_tilde
end



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

function norestart(opts...) end


"""
    perturbRestart!(x_tilde::Vector{Float64}, x_overline::Vector{Float64},
        indices::Vector{Int})

This perturbation replace the usual perturbation when a restart is needed, it is
the restart chosen in The feasibility pump 2005 by M. Fischetti, F. Glover and
A. Lodi
"""
function perturbRestart!(
        x_tilde::Vector{Float64},
        x_overline::Vector{Float64},
        indices::Vector{Int},
        opts...
    )
    nb_indices = length(indices)
    p = rand(nb_indices) .- 0.3

    for i in 1:nb_indices
        if abs(x_tilde[indices[i]] - x_overline[indices[i]]) + maximum([p[i], 0.0]) > 0.5
            x_tilde[indices[i]] = abs(x_tilde[indices[i]]-1)
        end
    end
    return x_tilde
end
