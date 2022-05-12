using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature
using Unitful

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, AbstractSetUp, AbstractWindDown,
AbstractAbiotic, _run_rule!, getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
_symmetric_grid, move!, get_neighbours, convert_coords
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
    #eco.spplist.species.traits.var[spp] = maximum(eco.spplist.species.traits.var)
    avgain = uconvert(NoUnits, prob * timestep) / size(eco.abundances.matrix, 2)
    gains = rand(rng, Poisson(avgain))
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
    prob::Float64
    length::Unitful.Time
    time::Unitful.Time
end

function getsoilwater(abenv::A) where A <: AbstractAbiotic
    return _getsoilwater(abenv.habitat)
end
function _getsoilwater(hab::ContinuousHab)
    return hab.matrix
end
function _getsoilwater(hab::HabitatCollection2)
    if eltype(hab.h1) <: Unitful.Volume
        return hab.h1.matrix
    elseif eltype(hab.h2) <: Unitful.Volume
        return hab.h2.matrix
    else
        error("No soil water volume layer")
    end
end

function Dry(location::Int64, prob::Float64, length::Unitful.Time)
    return Dry(location, prob, length, 1month)
end

mutable struct Rewet <: AbstractWindDown
    location::Int64
    prob::Float64
    length::Unitful.Time
    time::Unitful.Time
    max::Unitful.Volume
end
function Rewet(location::Int64, prob::Float64, length::Unitful.Time, max::Unitful.Volume)
    return Rewet(location, prob, length, 1month, max)
end

mutable struct WaterUse <: AbstractWindDown
    species::Int64
    location::Int64
    soil_moisture_frac::Float64
end

mutable struct LateralFlow <: AbstractSetUp
    location::Int64
    neighbours::Matrix{Int64}
    boundaries::Vector{Float64}
end
function LateralFlow(abenv::GridAbioticEnv, location::Int64)
    width = size(abenv.habitat, 1)
    x, y = convert_coords(location, width)
    neighbours = get_neighbours(getsoilwater(abenv), x, y, 8)
    boundary_length = sqrt.((neighbours[:, 1] .- x).^2 + (neighbours[:, 2] .- y).^2)
    return LateralFlow(location, neighbours, boundary_length)
end

mutable struct Ditch <: AbstractSetUp
    location::Int64
    prob::Float64
end

function _run_rule!(eco::Ecosystem, rule::LateralFlow, timestep::Unitful.Time)
    loc = rule.location
    neighbours = rule.neighbours
    boundaries = rule.boundaries
    hab = getsoilwater(eco.abenv)
    for i in 1:size(neighbours, 1)
        width = size(hab, 1)
        nei_loc = convert_coords(neighbours[i, 1], neighbours[i, 2], width)
        drainage = (hab[loc] - hab[nei_loc])/(20*boundaries[i])
        eco.cache.watermigration[loc] -= drainage
        eco.cache.watermigration[nei_loc] += drainage
    end
end

function _run_rule!(eco::Ecosystem, rule::WaterUse, timestep::Unitful.Time)
    spp = rule.species
    loc = rule.location
    hab = getsoilwater(eco.abenv)
    soil_moisture_frac = rule.soil_moisture_frac
    water_use = eco.spplist.species.requirement.energy[spp]
    area = getgridsize(eco)^2
    hab[loc] -= (soil_moisture_frac * water_use * area * eco.abundances.matrix[spp, loc])
end

function _run_rule!(eco::Ecosystem, rule::WaterFlux, timestep::Unitful.Time)
    loc = rule.location
    hab = getsoilwater(eco.abenv)
    area = getgridsize(eco)^2
    rainfall = eco.abenv.budget.matrix[loc] * area
    drainage = rule.prob * timestep * hab[loc]
    hab[loc] += rainfall - drainage
end

function _run_rule!(eco::Ecosystem, rule::Dry, timestep::Unitful.Time)
    loc = rule.location
    hab = getsoilwater(eco.abenv)
    drying = min(1, rule.time / rule.length) * rule.prob * hab[loc]
    hab[loc] -= drying
    rule.time += timestep
end

function _run_rule!(eco::Ecosystem, rule::Rewet, timestep::Unitful.Time)
    loc = rule.location
    hab = getsoilwater(eco.abenv)
    rewet = min(1, rule.time / rule.length) * rule.prob * (rule.max - hab[loc])
    hab[loc] += rewet
    rule.time += timestep
end

function changerate!(rule::Dry, newrate::DayType)
    return rule.prob = newrate
end