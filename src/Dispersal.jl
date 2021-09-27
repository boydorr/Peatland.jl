using EcoSISTEM: AbstractAbiotic
using EcoSISTEM
using SparseArrays

import EcoSISTEM: AbstractLookup


mutable struct PeatLookup <: AbstractLookup
    wind::SparseMatrixCSC{Float64, Int32}
    animal::SparseMatrixCSC{Float64, Int32}
    function PeatLookup(wind::SparseMatrixCSC{Float64, Int32}, animal::SparseMatrixCSC{Float64, Int32})
        all(0 .<= wind.nzval .<= 1) || error("Wind lookup values must be between 0 and 1")
        all(0 .<= animal.nzval .<= 1) || error("Animal lookup values must be between 0 and 1")
        return new(wind, animal)
    end
    function PeatLookup(Is::Vector{Int64}, Js::Vector{Int64}, total_size::Int64)
        Vs = fill(0.0, length(Is))
        wind = sparse(Int32.(Is), Int32.(Js), Vs, total_size, total_size)
        animal = sparse(Int32.(Is), Int32.(Js), Vs, total_size, total_size)
        return new(wind, animal)
    end
end

"""
  peatmove!(eco::AbstractEcosystem, id::Int64, pos::Int64, grd::Array{Int64, 2}, newvirus::Int64)

Function to calculate the movement of force of infection `id` from a given position in the landscape `pos`, using the lookup table found in the Ecosystem and updating the movement patterns on a cached grid, `grd`. The number of new virus is provided, so that movement only takes place as part of the generation process.
"""
function peatmove!(eco::AbstractEcosystem, id::Int64, pos::Int64, grd::Array{Float64, 2}, newseed::Int64)
    # Add in home movements
    wind = eco.lookup.wind
    if newseed > zero(newseed)
        for nzi in wind.colptr[pos]:(wind.colptr[pos+1]-1)
            grd[id, wind.rowval[nzi]] += newseed * wind.nzval[nzi]
        end
    end

    # Add in work movements
    animal = eco.lookup.animal
    if newseed > zero(newseed)
        for nzi in animal.colptr[pos]:(animal.colptr[pos+1]-1)
            grd[id, animal.rowval[nzi]] += newseed * animal.nzval[nzi]
        end
    end

    return eco
end


function _run_rule!(eco::Ecosystem, rule::WindDispersal)
    loc = getlocation(rule)
    spp = getspecies(rule)
    mov = eco.spplist.species.movement.kernels[spp]
    gridsize = getgridsize(eco)
    p = mov.thresh
    μ = (rule.plantheight) * (rule.windspeed / rule.terminalvelocity)
    σ_w = 0.5 * rule.windspeed
    σ = sqrt(rule.scaling * rule.canopyheight * 2(σ_w/rule.windspeed))
    λ = (rule.plantheight / σ)^2

    total_size = (size(eco.abenv.active, 1) * size(eco.abenv.active, 2))
    # Generate grid ids and x,y coords for active cells only
    grid_locs = 1:total_size
    activity = eco.abenv.active[1:end]
    grid_locs = grid_locs[activity]
    xys = convert_coords.(grid_locs, size(eco.abenv.active, 2))

    # Calculate lookup probabilities for each grid location
    locs, probs = genlookups(loc, grid_locs, xys, gridsize, p, λ, μ)
    eco.lookup.wind.nzval[locs] .= probs
end

function genlookups(from::Int64, to::Vector{Int64}, xys::Array{Tuple{Int64,Int64},1}, relsize::Unitful.Length, thresh::Float64, λ::Unitful.Length, μ::Unitful.Length)
    x, y = xys[to .== from][1]
    maxX = ceil(Int64, x + 1.0km/relsize) |> NoUnits; minX = ceil(Int64, x - 1.0km/relsize) |> NoUnits
    maxY = floor(Int64, y + 1.0km/relsize) |> NoUnits; minY = floor(Int64, y - 1.0km/relsize) |> NoUnits
    keep = [(i[1] <= maxX) & (i[2] <= maxY) & (i[1] >= minX) & (i[2] >= minY) for i in xys]
    to = to[keep]
    probs = [_lookup((x = x, y = y), (x = i[1], y = i[2]), relsize, _wald, λ, μ) for i in xys[keep]]
    keep = probs .> thresh
    probs = probs[keep]
    probs ./= sum(probs)
    return to[keep], probs
end

function _lookup(from::NamedTuple, to::NamedTuple, relSquareSize::Unitful.Length, dispersalfn::F, λ::Unitful.Length, μ::Unitful.Length) where {F<:Function}
    λ1 = λ / 1.0km |> NoUnits; 
    μ1 = μ / 1.0km |> NoUnits; 
    (relSquareSize /= 1.0km) |> NoUnits
    return hcubature(r -> dispersalfn(r, λ1, μ1),
      [from.y *relSquareSize - relSquareSize, from.x * relSquareSize - relSquareSize, to.y * relSquareSize - relSquareSize, to.x * relSquareSize - relSquareSize],
      [from.y * relSquareSize, from.x * relSquareSize, to.y * relSquareSize, to.x * relSquareSize],
      maxevals= 100, rtol = 0.01)[1]
end
