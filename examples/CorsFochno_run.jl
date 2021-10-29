using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots

file = "data/CorsFochno.tif"
cf = readfile(file, 261659.3m, 265409.3m, 289536.7m, 292756.7m)
active = Array(.!isnan.(Array(cf)))
heatmap(active)

grid = size(cf); individuals=5_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid); totalK = 0.1m^3/m^2
numMoss = 5; numShrub = 1; numSpecies = numMoss + numShrub
#Set up how much energy each species consumes
req1 = rand(Normal(1.0, 0.1), numMoss) .* m^3
req2 = rand(Normal(10.0, 0.1), numShrub) .* m^3
energy_vec = VolWaterRequirement([req1; req2])

# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(100.0m, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, NoBoundary())

# Create species list, including their temperature preferences, seed abundance and native status
pref1 = rand(Normal(100.0, 0.1), numMoss) .* m^3 
pref2 = rand(Normal(50.0, 0.1), numShrub) .* m^3 
opts = [pref1; pref2]
vars = fill(10.0m^3, numSpecies)
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
transitions = create_transition_list()
addtransition!(transitions, UpdateEnergy(update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
active_squares = findall(active[1:end])
for loc in active_squares
    for spp in eachindex(sppl.species.names) 
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, SeedDisperse(spp, loc))
    end
    addtransition!(transitions, Invasive(1, loc, 100.0/28days))
    addtransition!(transitions, Invasive(2, loc, 100.0/28days))
end

# Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)
eco.abundances.matrix[2, :] .= 0

# Run simulation
# Simulation Parameters
burnin = 5year; times1 = 5year; times2 = 5year; timestep = 1month;
record_interval = 1month; repeats = 1
lensim1 = length(0years:record_interval:times1)
lensim2 = length(0years:record_interval:times2)
abuns1 = generate_storage(eco, lensim1, 1)
abuns2 = generate_storage(eco, lensim2, 1)
abuns3 = generate_storage(eco, lensim2, 1)
# Burnin
@time simulate!(eco, burnin, timestep)

eco.abenv.habitat.change.rate = -1.0m^3/month

@time simulate_record!(abuns1, eco, times1, record_interval, timestep);
eco.abenv.habitat.change.rate = 0.0m^3/month
@time simulate_record!(abuns2, eco, times2, record_interval, timestep);
eco.abenv.habitat.change.rate = 1.0m^3/month
@time simulate_record!(abuns3, eco, times2, record_interval, timestep);

abuns = cat(abuns1, abuns2, abuns3, dims = 3)
abuns = reshape(abuns, (numSpecies, grid[1], grid[2], lensim1 + lensim2 *2, repeats))

plot(sum(abuns[1:numMoss, :, :, :, 1], dims = (1,2,3))[1, 1, 1, :], grid = false, label = "Moss")
plot!(sum(abuns[numMoss+1:end, :, :, :, 1], dims = (1,2,3))[1, 1, 1, :], label = "Shrub")
vline!([lensim1 + lensim2], color = :black, linestyle = :dash, label = "Intervention")
