using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature
using Unitful

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, AbstractSetUp, AbstractWindDown,
AbstractAbiotic, _run_rule!, getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
_symmetric_grid, move!, get_neighbours, convert_coords
const WaterTimeType = typeof(1.0m^3/day)
const VolType = typeof(1.0m^3)
const LengthType = typeof(1.0m)

abstract type AbstractPeatState <: AbstractStateTransition end
abstract type AbstractPeatPlace <: AbstractPlaceTransition end
abstract type AbstractPeatSetUp <: AbstractSetUp end
abstract type AbstractPeatWindDown <: AbstractWindDown end

"""
    Invasive <: AbstractStateTransition
Rule where an invasive species is added at a location with probability `prob`.
"""
mutable struct Invasive <: AbstractPeatState
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

"""
    WindDispersal <: AbstractPlaceTransition
Rule where a species is wind dispersed with an updating WALD kernel.
"""
mutable struct WindDispersal <: AbstractPeatPlace
    species::Int64
    location::Int64
    canopyheight::LengthType
    scaling::Float64
    plantheight::LengthType
    windspeed::Unitful.Velocity
    terminalvelocity::Unitful.Velocity
end

function _wald(r, λ, μ)
    d = ((r[3]-r[1])^2+(r[4]-r[2])^2)
    return (λ/(2*π*d^3))^0.5 * exp(-(λ*(d-μ)^2)/(2*μ^2*d))
end

"""
    WaterFlux <: AbstractSetUp
Rule where a particular location receives rainfall up to a maximum volume of `maxvol` and drains at a probability `prob`.
"""
mutable struct WaterFlux <: AbstractPeatSetUp
    location::Int64
    prob::DayType
    maxvol::VolType
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterFlux, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    maxvol = rule.maxvol
    area = getgridsize(eco)^2
    bud = eco.abenv.budget.matrix[loc]
    hab = eco.abenv.habitat.h1.matrix[loc]
    rainfall = bud * area
    drainage = rule.prob * timestep * hab
    hab = max(zero(typeof(drainage)), hab + rainfall - drainage)
    runoff = max(zero(typeof(maxvol)),  hab - maxvol)
    hab -= runoff
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterFlux, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    maxvol = rule.maxvol
    area = getgridsize(eco)^2
    rainfall = eco.abenv.budget.matrix[loc] * area
    drainage = rule.prob * timestep * eco.abenv.habitat.matrix[loc]
    eco.abenv.habitat.matrix[loc] = max(0.0m^3, eco.abenv.habitat.h1.matrix[loc] + rainfall - drainage)
    runoff = max(zero(typeof(maxvol)), eco.abenv.habitat.matrix[loc] - maxvol)
    eco.abenv.habitat.matrix[loc] -= runoff
end

"""
    Dry <: AbstractWindDown
Rule where a particular location dries out over a set length of time.
"""
mutable struct Dry <: AbstractPeatSetUp
    location::Int64
    prob::Float64
    length::Unitful.Time
    time::Unitful.Time
end

function Dry(location::Int64, prob::Float64, length::Unitful.Time)
    return Dry(location, prob, length, 1month)
end

function getprob(rule::R) where R <:AbstractSetUp
    return rule.prob
end
function getprob(rule::R) where R <:AbstractWindDown
    return rule.prob
end
function getlocation(rule::R) where R <:AbstractSetUp
    return rule.location
end
function getlocation(rule::R) where R <:AbstractWindDown
    return rule.location
end

function changerate!(rule::Dry, newrate::DayType)
    return rule.prob = newrate
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::Dry, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    drying = min(1, rule.time / rule.length) * rule.prob * eco.abenv.habitat.h1.matrix[loc]
    eco.abenv.habitat.h1.matrix[loc] -= drying
    rule.time += timestep
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::Dry, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    drying = min(1, rule.time / rule.length) * rule.prob * eco.abenv.habitat.matrix[loc]
    eco.abenv.habitat.matrix[loc] -= drying
    rule.time += timestep
end

"""
    Rewet <: AbstractWindDown
Rule where a particular location is rewet over a set length of time.
"""
mutable struct Rewet <: AbstractPeatSetUp
    location::Int64
    prob::Float64
    length::Unitful.Time
    time::Unitful.Time
    max::VolType
end
function Rewet(location::Int64, prob::Float64, length::Unitful.Time, max::Unitful.Volume)
    return Rewet(location, prob, length, 1month, max)
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::Rewet, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    rewet = min(1, rule.time / rule.length) * rule.prob * (rule.max - eco.abenv.habitat.h1.matrix[loc])
    eco.abenv.habitat.h1.matrix[loc] += rewet
    rule.time += timestep
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::Rewet, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    rewet = min(1, rule.time / rule.length) * rule.prob * (rule.max - eco.abenv.habitat.matrix[loc])
    eco.abenv.habitat.matrix[loc] += rewet
    rule.time += timestep
end

"""
    WaterUse <: AbstractStateTransition
Rule where a species at a location uses up the available soil moisture, given by `soil_moisture_frac`, according to its water use need.
"""
mutable struct WaterUse <: AbstractPeatState
    species::Int64
    location::Int64
    soil_moisture_frac::Float64
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterUse, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    spp = rule.species
    loc = rule.location
    soil_moisture_frac = rule.soil_moisture_frac
    hab = eco.abenv.habitat.h1.matrix[loc]
    abun = eco.abundances.matrix[spp, loc]
    water_needs = eco.spplist.species.requirement.energy[spp]
    area = getgridsize(eco)^2
    water_use = soil_moisture_frac * water_needs * area * abun
    hab = max(zero(typeof(water_use)), hab - water_use) 
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterUse, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    spp = rule.species
    loc = rule.location
    soil_moisture_frac = rule.soil_moisture_frac
    water_needs = eco.spplist.species.requirement.energy[spp]
    area = getgridsize(eco)^2
    water_use = soil_moisture_frac * water_needs * area * eco.abundances.matrix[spp, loc]
    eco.abenv.habitat.matrix[loc] = max(zero(typeof(water_use)), eco.abenv.habitat.matrix[loc] - water_use) 
end

"""
    LateralFlow <: AbstractWindDown
Rule where water flows from a grid square to its neighbours depending on the elevation and length of shared boundaries.
"""
mutable struct LateralFlow <: AbstractPeatSetUp
    location::Int64
    neighbours::Matrix{Int64}
    boundaries::Vector{Float64}
    elevation::Vector{Float64}
    prob::DayType
end

function LateralFlow(abenv::GridAbioticEnv, ele::ContinuousHab, location::Int64, rate::typeof(1.0/month))
    width = size(abenv.habitat, 1)
    x, y = convert_coords(location, width)
    neighbours = get_neighbours(abenv.habitat.h1.matrix, x, y, 8)
    boundary_length = sqrt.((neighbours[:, 1] .- x).^2 + (neighbours[:, 2] .- y).^2)
    elevation_diff = [(ele.matrix[neighbours[n, 1], neighbours[n, 2]] .- ele.matrix[x, y]) for n in 1:size(neighbours, 1)]
    return LateralFlow(location, neighbours, boundary_length, elevation_diff, rate)
end

function LateralFlow(abenv::GridAbioticEnv, location::Int64, rate::typeof(1.0/month))
    width = size(abenv.habitat, 1)
    x, y = convert_coords(location, width)
    neighbours = get_neighbours(abenv.habitat.matrix, x, y, 8)
    boundary_length = sqrt.((neighbours[:, 1] .- x).^2 + (neighbours[:, 2] .- y).^2)
    elevation_diff = fill(1, size(neighbours, 1))
    return LateralFlow(location, neighbours, boundary_length, elevation_diff, rate)
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::LateralFlow, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    neighbours = rule.neighbours
    boundaries = rule.boundaries
    for i in 1:size(neighbours, 1)
        width = size(eco.abenv.habitat, 1)
        nei_loc = convert_coords(neighbours[i, 1], neighbours[i, 2], width)
        drainage = rule.prob * timestep * rule.elevation[i] * (eco.abenv.habitat.h1.matrix[loc] - eco.abenv.habitat.h1.matrix[nei_loc])/(boundaries[i])
        eco.cache.watermigration[loc] -= drainage
        eco.cache.watermigration[nei_loc] += drainage
    end
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::LateralFlow, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    neighbours = rule.neighbours
    boundaries = rule.boundaries
    for i in 1:size(neighbours, 1)
        width = size(eco.abenv.habitat, 1)
        nei_loc = convert_coords(neighbours[i, 1], neighbours[i, 2], width)
        drainage = rule.prob * timestep * rule.elevation[i] * (eco.abenv.habitat.matrix[loc] - eco.abenv.habitat.matrix[nei_loc])/(boundaries[i])
        eco.cache.watermigration[loc] -= drainage
        eco.cache.watermigration[nei_loc] += drainage
    end
end