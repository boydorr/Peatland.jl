using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots

numSpecies = 2; grid = (10, 10); individuals=1000; area = 100.0*m^2; totalK = 10.0m^3/m^2

#Set up how much energy each species consumes
req1 = 1.0m^3; req2 = 10.0m^3;
energy_vec = VolWaterRequirement([req1, req2])

#Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, Torus())

# Create species list, including their temperature preferences, seed abundance and native status
pref1 = 100.0m^3; pref2 = 50.0m^3
opts = [pref1, pref2]
vars = [10.0m^3, 10.0m^3]
traits = GaussTrait(opts, vars)
native = [true, false]
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
abenv = peatlandAE(100.0m^3, grid, totalK, area)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0m^3)}()

# Create new transition list
transitions = TransitionList()
addtransition!(transitions, UpdateEnergy(update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for loc in eachindex(abenv.habitat.matrix)
    for spp in eachindex(sppl.species.names) 
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, SeedDisperse(spp, loc))
    end
    addtransition!(transitions, Invasive(2, loc, 1.0/month))
end

# Create ecosystem
eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
eco.abundances.matrix[2, :] .= 0

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
abuns = reshape(abuns, (numSpecies, grid[1], grid[2], lensim * 2, repeats))

# anim = @animate for i in 1:lensim*2
#    heatmap(abuns[1, :, :, i, 1], layout = (@layout [a b]), title = "Moss", clim = (0, 10))
#    heatmap!(abuns[2, :, :, i, 1], title = "Shrub", subplot = 2, clim = (0, 10))
# end
# gif(anim, fps = 30)

plot(sum(abuns[1, :, :, :, 1], dims = (1,2))[1, 1, :], grid = false, label = "Moss")
plot!(sum(abuns[2, :, :, :, 1], dims = (1,2))[1, 1, :], label = "Shrub")
vline!([lensim], color = :black, linestyle = :dash, label = "Intervention")
Plots.pdf("plots/Peatland_example.pdf")