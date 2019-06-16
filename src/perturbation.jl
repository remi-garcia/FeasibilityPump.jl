"""
Pertubation switchs binary variables: a switch means that a variable at value 0
will be switched at value 1 and a variable at value 1 will be switched at
value 0.
"""

function noperturb(opts...;compile::Bool=false) end


"""
    randomperturb!(xTilde::Vector{Float64}, _xOverline::Vector{Float64},
        indices::Vector{Int})

The random perturbation randomly switch variables. Each one has 50% chance to
be switched.
"""
function randomperturb!(
        xTilde::Vector{Float64},
        _xOverline::Vector{Float64},
        indices::Vector{Int},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [0.0]
    end

    for i in indices
        if bitrand()[1]
            xTilde[i] = abs(xTilde[i]-1)
        end
    end
    return xTilde
end


"""
    perturb!(xTilde::Vector{Float64}, xOverline::Vector{Float64},
        indices::Vector{Int}, TTmin::Int, TTmax::Int)

The perturbation switch the TT most fractional variables. TT is randomly picked
between TTmin to TTmax.
"""
function perturb!(
        xTilde::Vector{Float64},
        xOverline::Vector{Float64},
        indices::Vector{Int},
        TTmin::Int,
        TTmax::Int,
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [0.0]
    end

    nInd = length(indices)
    TT = rand(TTmin:TTmax)

    xTri = Vector{Tuple{Int, Float64}}()
    for i in indices
        push!(xTri, (i, abs(xOverline[i] - xTilde[i])))
    end
    sort!(xTri, by = x -> x[2], rev=true)

    for i in xTri[1:TT]
        xTilde[i[1]] = abs(xTilde[i[1]]-1)
    end
    return xTilde
end


"""
    frequencyperturb!(xTilde::Vector{Float64}, xOverline::Vector{Float64},
        indices::Vector{Int}, TTmin::Int, TTmax::Int)

The frequency perturbation switch TT variables. The variables are chosen by
frequency and TT is randomly picked between TTmin to TTmax..
"""
function frequencyperturb!(
        xTilde::Vector{Float64},
        xOverline::Vector{Float64},
        indices::Vector{Int},
        TTmin::Int,
        TTmax::Int,
        freqRound::Vector{Int},
        nIter::Int
        ;
        compile::Bool = false
    )
    if compile
        return [0.0]
    end

    nInd = length(indices)
    TT = rand(TTmin:TTmax)

    xTri = Vector{Tuple{Int, Float64}}()
    for i in indices
        if xTilde[i] < 0.5
            push!(xTri, (i, nIter - freqRound[i]))
        else
            push!(xTri, (i, freqRound[i]))
        end

    end
    sort!(xTri, by = x -> x[2], rev=true)

    for i in xTri[1:TT]
        xTilde[i[1]] = abs(xTilde[i[1]]-1)
    end
    return xTilde
end



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

function norestart(opts...;compile::Bool=false) end


"""
    perturbRestart!(xTilde::Vector{Float64}, xOverline::Vector{Float64},
        indices::Vector{Int})

This perturbation replace the usual perturbation when a restart is needed, it is
the restart chosen in The feasibility pump 2005 by M. Fischetti, F. Glover and
A. Lodi
"""
function perturbRestart!(
        xTilde::Vector{Float64},
        xOverline::Vector{Float64},
        indices::Vector{Int},
        opts...
        ;
        compile::Bool = false
    )
    if compile
        return [0.0]
    end

    nInd = length(indices)
    p = rand(nInd) .- 0.3

    for i in 1:nInd
        if abs(xTilde[indices[i]] - xOverline[indices[i]]) + maximum([p[i], 0.0]) > 0.5
            xTilde[indices[i]] = abs(xTilde[indices[i]]-1)
        end
    end
    return xTilde
end
