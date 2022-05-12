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

JLD2.@load("data/Peat_30_spp.jld2", peat_spp)

file = "data/CF_elevation.tif"
cf = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
active = Array(.!isnan.(Array(cf)))
heatmap(active)

grid = size(cf); individuals=5_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid); totalK = 0.01mm/m^2
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
pref1 = rand(Normal(100.0, 0.1), numMoss) .* m^3 
pref2 = rand(Normal(50.0, 0.1), numShrub) .* m^3 
opts = [pref1; pref2]
vars = fill(5.0m^3, numSpecies)
water_traits = GaussTrait(opts, vars)
ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
traits = TraitCollection2(water_traits, ele_traits)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
abenv = peatlandAE(cf, 100.0m^3, totalK) 

# Set relationship between species and environment (gaussian)
rel = additiveTR2(Gauss{typeof(1.0m^3)}(), Gauss{Float64}())

# Create new transition list
transitions = create_transition_list()
addtransition!(transitions, UpdateEnergy(update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
active_squares = findall(active[1:end])
for loc in active_squares
    for spp in eachindex(sppl.species.names) 
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, SeedDisperse(spp, loc))#, mean(plant_height), 0.3, plant_height[spp], 1.0m/s, 1.0m/s))
        addtransition!(transitions, WaterUse(spp, loc, 0.01))
        if spp > numMoss
            addtransition!(transitions, Invasive(spp, loc, 10.0/28days))
        end
    end
    addtransition!(transitions, LateralFlow(abenv, loc))
    addtransition!(transitions, WaterFlux(loc, 0.625/month))
end
# Add ditch drainage
file = "data/Ditches.tif"
ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
ditch_locations = findall(.!isnan.(ditches))
for d in ditch_locations
    loc = convert_coords(d[2], d[1], size(cf, 1))
    addtransition!(transitions, Dry(loc, 1.0, 1month))
end
# Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions, peatcache = true)

# Run simulation
# Simulation Parameters
burnin = 10year; times1 = 5year; times2 = 5year; timestep = 1month;
record_interval = 1month; repeats = 1
lensim = length(0years:record_interval:burnin)
lensim1 = length(0years:record_interval:times1)
lensim2 = length(0years:record_interval:times2)
abuns1 = generate_storage(eco, lensim1, 1)[:, :, :, 1]
abuns2 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
#abuns3 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
# Burnin
abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
@time simulate_record!(abuns, eco, burnin, record_interval, timestep)

sumabuns = Float64.(reshape(sum(abuns[:, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap(sumabuns')