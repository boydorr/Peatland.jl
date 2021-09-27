using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, _run_rule!,
getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
_symmetric_grid, move!

mutable struct Invasive <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::DayType
    function Invasive(species::Int64, location::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, prob)
    end
end

function _run_rule!(eco::Ecosystem, rule::Invasive, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    prob = getprob(rule)
    eco.spplist.species.traits.var[spp] = maximum(eco.spplist.species.traits.var)
    invasive_abun = eco.spplist.species.abun[spp]
    avgain = uconvert(NoUnits, prob * timestep) / size(eco.abundances.matrix, 2)
    gains = rand(Poisson(avgain))
    eco.abundances.matrix[spp, loc] += gains
end

mutable struct WindDispersal <: AbstractPlaceTransition
    species::Int64
    location::Int64
    canopyheight::Unitful.Length
    scaling::Float64
    plantheight::Unitful.Length
    windspeed::Unitful.Velocity
    terminalvelocity::Unitful.Velocity
end

# function _run_rule!(eco::Ecosystem, rule::WindDispersal)
#     loc = getlocation(rule)
#     spp = getspecies(rule)
#     hab = eco.abenv.habitat
#     mov = eco.spplist.species.movement.kernels[spp]
#     gridsize = getgridsize(eco)
#     m = maximum(_getdimension(hab))
#     p = mov.thresh
#     μ = (rule.plantheight / gridsize) * (rule.windspeed / rule.terminalvelocity)
#     σ_w = 0.5 * rule.windspeed
#     σ = sqrt(rule.scaling * rule.canopyheight * 2(σ_w/rule.windspeed))
#     λ = (rule.plantheight / σ)^2 / gridsize
#     lookup_dat = _lookup(m, p, _wald, uconvert(NoUnits, λ), uconvert(NoUnits, μ))
#     eco.lookup.species[spp] = Lookup(lookup_dat)
#     if eco.abenv.active[loc]
#         move!(eco, eco.spplist.species.movement, loc, spp, eco.cache.netmigration, eco.cache.seedbank[spp, loc])
#     end
# end

function _wald(r, λ, μ)
    d = ((r[3]-r[1])^2+(r[4]-r[2])^2)
    return (λ/(2*π*d^3))^0.5 * exp(-(λ*(d-μ)^2)/(2*μ^2*d))
   end

# function _lookup(maxGridSize::Int64,
#                 pThresh::Float64, dispersalfn::F, λ::Float64, μ::Float64) where {F<:Function}
#   # Create empty array
#   lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])
#   # Loop through directions until probability is below threshold
#   k = 0
#   m = 0
#   count = 0
#   while (k <= maxGridSize && m <= maxGridSize)
#     count = count + 1
#     calc_prob = hcubature(r -> dispersalfn(r, λ, μ),
#       [0, 0, k, m],
#       [1, 1, (k+1), (m+1)],
#       maxevals= 10000)[1]
#     if m == 0 && calc_prob < pThresh
#       break
#     end
#     if count == 1
#       push!(lookup_tab, [k m calc_prob])
#       k = k + 1
#     elseif (calc_prob > pThresh && m <= k)
#       push!(lookup_tab, [k m calc_prob])
#       m = m + 1
#     else
#       m = 0
#       k = k + 1
#     end
#   end
#   # If no probabilities can be calculated, threshold is too high
#   nrow(lookup_tab) != 0 || error("probability threshold too high")
#   # Find all other directions
#   lookup_tab = _symmetric_grid(lookup_tab)
#   #info(sum(lookup_tab[:, 3]))
#   # Normalise
#   lookup_tab[!, :Prob] = lookup_tab[!, :Prob]/sum(lookup_tab[!, :Prob])
#   return lookup_tab
# end
