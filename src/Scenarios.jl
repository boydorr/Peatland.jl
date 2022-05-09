using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature
using Unitful

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, AbstractSetUp, AbstractWindDown,
_run_rule!, getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
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
    eco.spplist.species.traits.var[spp] = maximum(eco.spplist.species.traits.var)
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
    width = size(abenv.habitat.matrix, 1)
    x, y = convert_coords(location, width)
    neighbours = get_neighbours(abenv.habitat.matrix, x, y, 8)
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
    for i in 1:size(neighbours, 1)
        width = size(eco.abenv.habitat.matrix, 1)
        nei_loc = convert_coords(neighbours[i, 1], neighbours[i, 2], width)
        drainage = (eco.abenv.habitat.matrix[loc] - eco.abenv.habitat.matrix[nei_loc])/(20*boundaries[i])
        eco.cache.watermigration[loc] -= drainage
        eco.cache.watermigration[nei_loc] += drainage
    end
end

function _run_rule!(eco::Ecosystem, rule::WaterUse, timestep::Unitful.Time)
    spp = rule.species
    loc = rule.location
    soil_moisture_frac = rule.soil_moisture_frac
    water_use = eco.spplist.species.requirement.energy[spp]
    area = getgridsize(eco)^2
    eco.abenv.habitat.matrix[loc] -= (soil_moisture_frac * water_use * area * eco.abundances.matrix[spp, loc])
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
    drying = min(1, rule.time / rule.length) * rule.prob * eco.abenv.habitat.matrix[loc]
    eco.abenv.habitat.matrix[loc] -= drying
    rule.time += timestep
end

function _run_rule!(eco::Ecosystem, rule::Rewet, timestep::Unitful.Time)
    loc = rule.location
    rewet = min(1, rule.time / rule.length) * rule.prob * (rule.max - eco.abenv.habitat.matrix[loc])
    eco.abenv.habitat.matrix[loc] += rewet
    rule.time += timestep
end

function changerate!(rule::Dry, newrate::DayType)
    return rule.prob = newrate
end