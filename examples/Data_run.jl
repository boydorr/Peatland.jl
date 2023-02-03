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

file = "data/CorsFochno.tif"
cf = readfile(file, 261659.3m, 265409.3m, 289536.7m, 292756.7m)
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
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
abenv = peatlandAE(100.0m^3, grid, totalK, area, active)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0m^3)}()

# Create new transition list
transitions = TransitionList()
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
    κ = 0.01/month
    λ = 0.01/month

    addtransition!(transitions, LateralFlow(loc, κ, λ))
    addtransition!(transitions, WaterFlux(loc, 0.625/month))
end

# Create ecosystem
eco = PeatSystem(sppl, abenv, rel, transitions = transitions)



# envs = zeros(lensim)
# for i in 1:lensim
#     EcoSISTEM.update!(eco, timestep, transitions)
#     envs[i] = ustrip.(mean(eco.abenv.habitat.matrix[eco.abenv.active]))
# end
# plot(envs)

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
# plot(sum(abuns[1:numMoss, :, :, :, 1], dims = (1,2))[1, 1, :], grid = false, label = "",
# layout = (2, 1), title = "Moss")
# plot!(sum(abuns[numMoss+1:end, :, :, :, 1], dims = (1,2))[1, 1, :], label = "",
# subplot = 2, color = 2, grid = false, title = "Shrub")
# vline!([lensim1 + lensim2], color = :black, linestyle = :dash, label = "Intervention")
# vline!([lensim1 + lensim2], color = :black, linestyle = :dash, label = "",
# subplot = 2)

for loc in active_squares
    addtransition!(transitions, Peatland.Dry(loc, 0.5, 1year))
end

@time simulate_record!(abuns1, eco, times1, record_interval, timestep);

# rules = eco.transitions.winddown[typeof.(eco.transitions.winddown) .== Dry]
# changerate!.(rules, 0.0/30days)
for loc in active_squares
    addtransition!(transitions, Peatland.Rewet(loc, 0.5, 1year, 150.0m^3))
end

@time simulate_record!(abuns2, eco, times2, record_interval, timestep);
# for loc in active_squares
#     addtransition!(transitions, Peatland.Rewet(loc, 0.5/30days, 100.0m^3))
# end
# @time simulate_record!(abuns3, eco, times2, record_interval, timestep, save = true);

abuns = cat(abuns, abuns1, abuns2,  dims = 3)
abuns = reshape(abuns, (numSpecies, grid[1], grid[2], lensim + lensim1 + lensim2, repeats))

times = collect(1:231) ./ 12
p = plot(times, sum(abuns[1:numMoss, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :]./10, grid = false, label = "",
layout = (3, 1), title = "Moss", size = (1000, 1200))
p = plot!(times, sum(abuns[shrubs, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :], label = "",
subplot = 2, color = 2, grid = false, title = "Shrubs, Herbs and Grasses")
p = plot!(times, sum(abuns[trees, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :], label = "",
subplot = 3, color = 3, grid = false, title = "Trees")
for i in 1:3
    p = vline!(times, [(lensim-12)/12], color = :black, linestyle = :dot, label = "Drying", subplot = i)
    p = vline!(times, [(lensim+lensim1-12)/12], color = :black, linestyle = :dash, label = "Intervention", subplot = i)
end
mean1 = mean(sum(abuns[mosses, :, :, 13:lensim, 1], dims = (1,2,3))[1, 1, 1, :]./10)
mean2 = mean(sum(abuns[mosses, :, :, (lensim+lensim1+13):end, 1], dims = (1,2,3))[1, 1, 1, :]./10)
mean3 = mean(sum(abuns[shrubs, :, :, 13:lensim, 1], dims = (1,2,3))[1, 1, 1, :])
mean4 = mean(sum(abuns[shrubs, :, :, (lensim+lensim1+13):end, 1], dims = (1,2,3))[1, 1, 1, :])
mean5 = mean(sum(abuns[trees, :, :, 13:lensim, 1], dims = (1,2,3))[1, 1, 1, :])
mean6 = mean(sum(abuns[trees, :, :, (lensim+lensim1+13):end, 1], dims = (1,2,3))[1, 1, 1, :])
p = plot!(collect(1:lensim)/12, fill(mean1, 121), color = 1, linestyle = :dot, label = "Start mean", subplot = 1)
p = plot!(collect((lensim+lensim1):243)/12, fill(mean2, 62), color = 1, linestyle = :dot, label = "Final mean", subplot = 1)
p = plot!(collect(1:lensim)/12, fill(mean3, 121), color = 2, linestyle = :dot, label = "Start mean", subplot = 2)
p = plot!(collect((lensim+lensim1):243)/12, fill(mean4, 62), color = 2, linestyle = :dot, label = "Final mean", subplot = 2)
p = plot!(collect(1:lensim)/12, fill(mean5, 121), color = 3, linestyle = :dot, label = "Start mean", subplot = 3)
p = plot!(collect((lensim+lensim1):243)/12, fill(mean6, 62), color = 3, linestyle = :dot, label = "Final mean", subplot = 3)
display(p)
Plots.pdf("plots/Water_cycle2.pdf")
Plots.png("plots/Water_cycle2.png")