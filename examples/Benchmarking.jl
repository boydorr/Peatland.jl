using Pkg
Pkg.activate("..")
Pkg.instantiate()
using BenchmarkTools
using Profile
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots
using JLD2
using DataFrames
using Missings
using CSV
using Distances
using Shapefile
using AxisArrays

JLD2.@load("../data/Peat_30_spp.jld2", peat_spp)

file = "../data/CF_elevation.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
active = Array(.!isnan.(Array(cf)))

grid = (40, 50); individuals=10_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid);
numMoss = 5; numShrub = 28; numSpecies = numMoss + numShrub

mosses = 1:numMoss
shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss

height = peat_spp[!, :Plant_height]
shrub_size = height .* 1.0m^2
moss_height = fill(1.0, numMoss)
moss_size = moss_height .* 1.0m^2
plant_height = [moss_height; height] .* m
#Set up how much energy each species consumes
req1 = moss_size .* rand(Normal(1.0, 0.01), numMoss) .* mm ./m^2
req2 = shrub_size .* rand(Normal(10.0, 0.1), numShrub) .* mm ./m^2
energy_vec = WaterRequirement([req1; req2])

# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(100.0m, 1e-3), numSpecies)
movement = BirthOnlyMovement(kernel, NoBoundary())

# Create species list, including their temperature preferences, seed abundance and native status
pref1 = rand(Normal(70.0, 0.1), numMoss) .* m^3 
pref2 = rand(Normal(30.0, 0.1), numShrub) .* m^3 
opts = [pref1; pref2]
vars = [fill(5.0m^3, numMoss); fill(10.0m^3, numShrub)]
water_traits = GaussTrait(opts, vars)
ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
soil_pref1 = fill([1], numMoss); soil_pref2 = fill([0,1], numShrub)
soil_pref = soiltrait([soil_pref1; soil_pref2])
traits = TraitCollection3(water_traits, ele_traits, soil_pref)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = [fill(div(individuals, numMoss), numMoss); fill(0, numShrub)]
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
ele = fill(10.0, grid)
ele = AxisArray(ele, Axis{:latitude}(1.0m:10.0m:400m), Axis{:longitude}(1.0m:10.0m:500m))

bud = WaterBudget(fill(100.0mm, grid))
soil = fill(1, grid)
abenv = peatlandAE(ele, soil, 100.0m^3, bud) 
# abenv.habitat.h1.matrix ./= abenv.habitat.h2.matrix
# centre_point = (180, 270)
# dists = [euclidean(centre_point, i) for i in Tuple.(findall(active))[:]]
# active_squares = findall(active)
# for i in eachindex(active_squares)
#     abenv.habitat.h1.matrix[active_squares[i]] += 3 .* exp.(-1 .* dists[i]/70) * abenv.habitat.h1.matrix[active_squares[i]]
# end 
# abenv.habitat.h1.matrix[isinf.(abenv.habitat.h1.matrix)] .= 100.0m^3
# abenv.habitat.h1.matrix[isnan.(abenv.habitat.h1.matrix)] .= 0.0m^3

# Set relationship between species and environment (gaussian)
rel = multiplicativeTR3(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}(), soilmatch{Int64}())

# Create new transition list
transitions = TransitionList()
addtransition!(transitions, UpdateEnergy(update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
active_squares = findall(active[1:end])
peat_squares = findall(soil .== 1)
for loc in eachindex(abenv.habitat.h1.matrix)
    for spp in eachindex(sppl.species.names) 
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, GaussDispersal(spp, loc))
        addtransition!(transitions, WaterUse(spp, loc, 0.01))
        if spp > numMoss
            addtransition!(transitions, Invasive(spp, loc, 10.0/28days))
        end
    end
    if loc ∈ peat_squares
        drainage = 0.00001/month; flow = 0.0001/month
    else
        drainage = 0.006/month; flow = 0.5/month
    end
    addtransition!(transitions, LateralFlow(abenv, abenv.habitat.h2, loc, flow))
    addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
end

eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
Profile.clear_malloc_data()
# Run simulation
# Simulation Parameters
burnin = 1month; times = 2month; timestep = 1month;
simulate!(eco, burnin, timestep)
Profile.clear_malloc_data()
simulate!(eco, times, timestep)