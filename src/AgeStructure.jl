using EcoSISTEM: getdestination
using EcoSISTEM
import EcoSISTEM: AbstractStateTransition, _run_rule!, TimeType


mutable struct Ageing <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::TimeType
end

function _run_rule!(eco::Ecosystem, rule::Ageing, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if (eco.abenv.active[loc])
        ageprob = getprob(rule) * timestep
        newageprob = 1.0 - exp(-ageprob)
        age_up = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newageprob))
        eco.abundances.matrix[spp, loc] -= age_up
        eco.abundances.matrix[dest, loc] += age_up
    end
end

