using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature
using Unitful

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, AbstractSetUp, AbstractWindDown,
_run_rule!, getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
_symmetric_grid, move!
const WaterTimeType = typeof(1.0m^3/day)

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

function _wald(r, λ, μ)
    d = ((r[3]-r[1])^2+(r[4]-r[2])^2)
    return (λ/(2*π*d^3))^0.5 * exp(-(λ*(d-μ)^2)/(2*μ^2*d))
   end

mutable struct WaterFlux <: AbstractSetUp
    location::Int64
    prob::DayType
end

mutable struct Dry <: AbstractWindDown
    location::Int64
    prob::DayType
end

mutable struct WaterUse <: AbstractWindDown
    species::Int64
    location::Int64
end

function _run_rule!(eco::Ecosystem, rule::WaterUse, timestep::Unitful.Time)
    spp = rule.species
    loc = rule.location
    water_use = eco.spplist.species.requirement.energy[spp]
    area = getgridsize(eco)^2
    eco.abenv.habitat.matrix[loc] -= (water_use * area * eco.abundances.matrix[spp, loc])
end

function _run_rule!(eco::Ecosystem, rule::WaterFlux, timestep::Unitful.Time)
    loc = rule.location
    area = getgridsize(eco)^2
    rainfall = eco.abenv.budget.matrix[loc] * area
    drainage = rule.prob * timestep * eco.abenv.habitat.matrix[loc]
    eco.abenv.habitat.matrix[loc] += rainfall - drainage
end

function _run_rule!(eco::Ecosystem, rule::Dry, timestep::Unitful.Time)
    loc = rule.location
    drying = rule.prob * timestep * eco.abenv.habitat.matrix[loc]
    eco.abenv.habitat.matrix[loc] -= drying
end

function changerate!(rule::Dry, newrate::DayType)
    return rule.prob = newrate
end