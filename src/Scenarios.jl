using EcoSISTEM.Units
using Distributions
using DataFrames
using HCubature
using Unitful

import EcoSISTEM: DayType, AbstractPlaceTransition, AbstractStateTransition, AbstractSetUp, AbstractWindDown,
AbstractAbiotic, _run_rule!, getspecies, getlocation, getprob, TimeType, Lookup, _getdimension,
_symmetric_grid, move!, get_neighbours, convert_coords
const WaterTimeType = typeof(1.0m^3/day)
const AreaTimeType = typeof(1.0m^2/day)
const VolType = typeof(1.0m^3)
const LengthType = typeof(1.0m)
const InvVolType = typeof(1/1.0m^3)

abstract type AbstractPeatState <: AbstractStateTransition end
abstract type AbstractPeatPlace <: AbstractPlaceTransition end
abstract type AbstractPeatSetUp <: AbstractSetUp end
abstract type AbstractPeatWindDown <: AbstractWindDown end

"""
    Invasive <: AbstractStateTransition
Rule where an invasive species is added at a location with probability `prob`.
"""
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

"""
    WindDispersal <: AbstractPlaceTransition
Rule where a species is wind dispersed with an updating WALD kernel.
"""
mutable struct WindDispersal <: AbstractPlaceTransition
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
mutable struct WaterFlux <: AbstractSetUp
    location::Int64
    fmax::DayType
    uptake_rate::InvVolType
    function WaterFlux(location::Int64, fmax::T, uptake_rate::InvVolType) where T
        fmaxnew = uconvert(unit(DayType), fmax)
        new(location, fmaxnew, uptake_rate)
    end
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterFlux, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    area = getgridsize(eco)^2
    bud = eco.abenv.budget.matrix[loc]
    hab = eco.abenv.habitat.h1.matrix[loc]
    rainfall = bud * area
    eco.cache.surfacewater[loc] += rainfall
    infiltration = rule.fmax * (1 - hab) * timestep * eco.cache.surfacewater[loc]
    eco.cache.surfacewater[loc] = max(zero(typeof(infiltration)), eco.cache.surfacewater[loc] - infiltration)
    eco.abenv.habitat.h1.matrix[loc] += infiltration * rule.uptake_rate
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterFlux, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    area = getgridsize(eco)^2
    bud = eco.abenv.budget.matrix[loc]
    hab = eco.abenv.habitat.matrix[loc]
    rainfall = bud * area
    eco.cache.surfacewater[loc] += rainfall
    infiltration = rule.fmax * (1 - hab) * timestep * eco.cache.surfacewater[loc]
    eco.cache.surfacewater[loc] = max(zero(typeof(infiltration)), eco.cache.surfacewater[loc] - infiltration)
    eco.abenv.habitat.matrix[loc] += infiltration * rule.uptake_rate
end

"""
    Dry <: AbstractWindDown
Rule where a particular location dries out over a set length of time.
"""
mutable struct Dry <: AbstractSetUp
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
mutable struct Rewet <: AbstractSetUp
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
mutable struct WaterUse <: AbstractStateTransition
    species::Int64
    location::Int64
    uptake_rate::InvVolType
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterUse, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    spp = rule.species
    loc = rule.location
    # uptake_rate = rule.uptake_rate
    hab = eco.abenv.habitat.h1.matrix[loc]
    abun = eco.abundances.matrix[spp, loc]
    water_needs = eco.spplist.species.requirement.energy[spp]
    bud = eco.abenv.budget.matrix[loc]
    uptake_rate = water_needs/bud
    # area = getgridsize(eco)^2
    water_use = uptake_rate * abun * hab
    eco.abenv.habitat.h1.matrix[loc] = max(zero(typeof(water_use)), hab - water_use) 
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::WaterUse, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    spp = rule.species
    loc = rule.location
    # uptake_rate = rule.uptake_rate
    hab = eco.abenv.habitat.matrix[loc]
    abun = eco.abundances.matrix[spp, loc]
    water_needs = eco.spplist.species.requirement.energy[spp]
    bud = eco.abenv.budget.matrix[loc]
    uptake_rate = water_needs/bud
    # area = getgridsize(eco)^2
    water_use = uptake_rate * abun * hab
    eco.abenv.habitat.matrix[loc] = max(zero(typeof(water_use)), hab - water_use) 
end

"""
    LateralFlow <: AbstractWindDown
Rule where water flows from a grid square to its neighbours depending on the elevation and length of shared boundaries.
"""
mutable struct LateralFlow <: AbstractStateTransition
    location::Int64
    κ::AreaTimeType
    λ::AreaTimeType
end

function getspecies(rule::LateralFlow)
    return 1
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::LateralFlow, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    loc = rule.location
    x, y = convert_coords(loc, size(eco.abenv.habitat, 1))
    gs = getgridsize(eco)
    maxX = size(eco.abenv.habitat, 1)
    maxY = size(eco.abenv.habitat, 2)
    if (x > 1) && (x < maxX) && (y > 1) && (y < maxY)
        update_ghostcells!(eco.abenv.habitat.h1.matrix)

        u1 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x + 1, y] - eco.abenv.habitat.h2.matrix[x - 1, y]) / 2gs
        u2 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x + 1, y] - 2 * eco.abenv.habitat.h2.matrix[x, y] + eco.abenv.habitat.h2.matrix[x - 1, y]) / gs^2
        v1 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x, y + 1] - eco.abenv.habitat.h2.matrix[x, y - 1]) / 2gs
        v2 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x, y + 1] - 2 * eco.abenv.habitat.h2.matrix[x, y] + eco.abenv.habitat.h2.matrix[x, y - 1]) / gs^2
        advection_x1 =  u1 * (eco.cache.surfacewater[x + 1, y] - eco.cache.surfacewater[x - 1, y]) / 2gs
        advection_x2 = u2 * eco.cache.surfacewater[x, y]
        advection_y1 = v1 * (eco.cache.surfacewater[x, y + 1] - eco.cache.surfacewater[x, y - 1]) / 2gs
        advection_y2 = v2 * eco.cache.surfacewater[x, y]
        advection = advection_x1 + advection_x2 + advection_y1 + advection_y2

        diffusion_x = (eco.cache.surfacewater[x + 1, y] - 2*eco.cache.surfacewater[x, y] + eco.cache.surfacewater[x - 1, y]) / gs^2
        diffusion_y = (eco.cache.surfacewater[x, y + 1] - 2*eco.cache.surfacewater[x, y] + eco.cache.surfacewater[x, y - 1]) / gs^2
        diffusion = rule.κ * timestep * (diffusion_x + diffusion_y)

        eco.cache.surfacemigration[loc] += diffusion + advection 
    end
end

# function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::LateralFlow, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
#     loc = rule.location
#     x, y = convert_coords(loc, size(eco.abenv.habitat, 1))
#     gs = getgridsize(eco)
#     maxX = size(eco.abenv.habitat, 1)
#     maxY = size(eco.abenv.habitat, 2)
#     if (x > 1) && (x < maxX) && (y > 1) && (y < maxY)
#         update_ghostcells!(eco.abenv.habitat.h1.matrix)

#         u1 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x + 1, y] - eco.abenv.habitat.h2.matrix[x - 1, y]) / 2gs
#         u2 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x + 1, y] - 2 * eco.abenv.habitat.h2.matrix[x, y] + eco.abenv.habitat.h2.matrix[x - 1, y]) / gs^2
#         v1 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x, y + 1] - eco.abenv.habitat.h2.matrix[x, y - 1]) / 2gs
#         v2 = rule.λ * timestep * (eco.abenv.habitat.h2.matrix[x, y + 1] - 2 * eco.abenv.habitat.h2.matrix[x, y] + eco.abenv.habitat.h2.matrix[x, y - 1]) / gs^2
#         advection_x1 =  u1 * (eco.abenv.habitat.h1.matrix[x + 1, y] - eco.abenv.habitat.h1.matrix[x - 1, y]) / 2gs
#         advection_x2 = u2 * eco.abenv.habitat.h1.matrix[x, y]
#         advection_y1 = v1 * (eco.abenv.habitat.h1.matrix[x, y + 1] - eco.abenv.habitat.h1.matrix[x, y - 1]) / 2gs
#         advection_y2 = v2 * eco.abenv.habitat.h1.matrix[x, y]
#         advection = advection_x1 + advection_x2 + advection_y1 + advection_y2

#         diffusion_x = (eco.abenv.habitat.h1.matrix[x + 1, y] - 2*eco.abenv.habitat.h1.matrix[x, y] + eco.abenv.habitat.h1.matrix[x - 1, y]) / gs^2
#         diffusion_y = (eco.abenv.habitat.h1.matrix[x, y + 1] - 2*eco.abenv.habitat.h1.matrix[x, y] + eco.abenv.habitat.h1.matrix[x, y - 1]) / gs^2
#         diffusion = rule.κ * timestep * (diffusion_x + diffusion_y)

#         eco.cache.watermigration[loc] += diffusion + advection
#     end
# end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::LateralFlow, timestep::Unitful.Time) where {A, B, H <: ContinuousHab}
    loc = rule.location
    x, y = convert_coords(loc, size(eco.abenv.habitat, 1))
    maxX = size(eco.abenv.habitat, 1)
    maxY = size(eco.abenv.habitat, 2)
    if (x > 1) && (x < maxX) && (y > 1) && (y < maxY)
        update_ghostcells!(eco.abenv.habitat.matrix)
        diffusion_x = (eco.cache.surfacewater[x + 1, y] - 2*eco.cache.surfacewater[x, y] + eco.cache.surfacewater[x - 1, y]) / gs^2
        diffusion_y = (eco.cache.surfacewater[x, y + 1] - 2*eco.cache.surfacewater[x, y] + eco.cache.surfacewater[x, y - 1]) / gs^2
        diffusion = rule.κ * timestep * (diffusion_x + diffusion_y)

        eco.cache.surfacemigration[loc] += diffusion
    end
end

mutable struct Drainage <: AbstractStateTransition
    location::Int64
    drainage::WaterTimeType
end

function _run_rule!(eco::Ecosystem{A, GridAbioticEnv{H, B}}, rule::Drainage, timestep::Unitful.Time) where {A, B, H <: Union{HabitatCollection2, HabitatCollection3}}
    hab = eco.cache.surfacewater[loc]
    drainage = rule.drainage * timestep * hab
    eco.cache.surfacewater[loc] = max(zero(typeof(drainage)), hab + rainfall - drainage)
end
