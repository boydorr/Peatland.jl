using EcoSISTEM.Units
using Distributions

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, _run_rule!,
getspecies, getlocation, getprob

mutable struct Invasive <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::DayType
end


function _run_rule!(eco::Ecosystem, rule::Invasive, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    prob = getprob(rule)
    eco.spplist.species.traits.var[spp] = maximum(eco.spplist.species.traits.var)
    invasive_abun = eco.spplist.species.abun[spp]
    avgain = uconvert(NoUnits, prob * timestep)
    gains = rand(Poisson(avgain))
    pos = 1:size(eco.abundances.matrix, 2)
    smp = sample(pos, gains)
    eco.abundances.grid[spp, smp] .+= 1
end
