using EcoSISTEM.Units
using Distributions

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, _run_rule!,
getspecies, getlocation, getprob, TimeType

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
