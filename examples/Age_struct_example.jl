using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots

numSpecies = 2; numAgeClasses = 2; grid = (10, 10); individuals=1000; area = 100.0*m^2; totalK = 10.0m^3/m^2

#Set up how much energy each species consumes
req1 = fill(1.0m^3, numAgeClasses); req2 = fill(10.0m^3, numAgeClasses);
energy_vec = VolWaterRequirement([req1; req2])

#Set rates for birth and death
birth = 0.5/year
death = 0.5/year
ageing = 1.0/month
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies * numAgeClasses)
movement = BirthOnlyMovement(kernel, Torus())

# Create species list, including their temperature preferences, seed abundance and native status
pref1 = fill(100.0m^3, numAgeClasses); pref2 = fill(50.0m^3, numAgeClasses)
opts = [pref1; pref2]
vars = fill(10.0m^3, numAgeClasses * numSpecies)
traits = GaussTrait(opts, vars)
native = [fill(true, numAgeClasses); fill(false, numAgeClasses)]
abun = fill(div(individuals, numSpecies * numAgeClasses), numSpecies * numAgeClasses)
sppl = SpeciesList(numSpecies * numAgeClasses, traits, abun, energy_vec,
    movement, param, native)
sppl.species.names = ["Moss_gametophyte", "Moss_sporophyte", "Shrub_seedling", "Shrub_mature"]

# Create abiotic environment - even grid of one temperature
abenv = peatlandAE(100.0m^3, grid, totalK, area)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0m^3)}()

# Age by species categories
cat_idx = reshape(1:(numSpecies * numAgeClasses), numAgeClasses, numSpecies)
identities = (young = cat_idx[1, :], mature = cat_idx[2, :])

sppl.params.death[identities.young] *= 2
# Create new transition list
transitions = create_transition_list()
addtransition!(transitions, UpdateEnergy(update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for loc in eachindex(abenv.habitat.matrix)
    # Death & ageing can happen to any category
    for spp in eachindex(sppl.species.names)
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
    end
    # Add invasive shrub seedlings
    addtransition!(transitions, Invasive(3, loc, 1.0/month))

    # Only young plants can grow
    for young_spp in identities.young
        addtransition!(transitions, Ageing(young_spp, loc, young_spp + 1, ageing))
    end
    # Only mature plants produce seed
    for mature_spp in identities.mature
        addtransition!(transitions, GenerateSeed(mature_spp, loc, mature_spp - 1, sppl.params.birth[mature_spp]))   
        addtransition!(transitions, SeedDisperse(mature_spp, loc))
    end
end

# Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)
eco.abundances.matrix[cat_idx[:, 2], :] .= 0

# Run simulation
# Simulation Parameters
burnin = 1year; times = 5year; timestep = 1day;
record_interval = 1day; repeats = 1
lensim = length(0years:record_interval:times)
abuns1 = generate_storage(eco, lensim, 1)
abuns2 = generate_storage(eco, lensim, 1)
# Burnin
@time simulate!(eco, burnin, timestep)
@time simulate_record!(abuns1, eco, times, record_interval, timestep)

eco.abenv.habitat.change.rate = -1.0m^3/month
@time simulate_record!(abuns2, eco, times, record_interval, timestep)

abuns = cat(abuns1, abuns2, dims = 3)
abuns = reshape(abuns, (numSpecies * numAgeClasses, grid[1], grid[2], lensim * 2, repeats))

# anim = @animate for i in 1:lensim*2
#    heatmap(abuns[1, :, :, i, 1], layout = (@layout [a b]), title = "Moss", clim = (0, 10))
#    heatmap!(abuns[2, :, :, i, 1], title = "Shrub", subplot = 2, clim = (0, 10))
# end
# gif(anim, fps = 30)

plot(grid = false)
for i in eachindex(sppl.species.names)
display(plot!(sum(abuns[i, :, :, :, 1], dims = (1,2))[1, 1, :], label = sppl.species.names[i]))
end
vline!([lensim], color = :black, linestyle = :dash, label = "Drying")

Plots.pdf("plots/Peatland_age_example.pdf")

@time simulate!(eco, 1year, 1day)