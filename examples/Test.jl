using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity
using DataPipeline
using Plots

handle = initialise()
# Set up initial parameters for ecosystem
numSpecies = 10; grid = (10, 10); req= 10.0kJ; individuals=100_000; area = 100.0*km^2; totalK = 1_000.0kJ/km^2

# Set up how much energy each species consumes
# energy_vec = SolarRequirement(collect(1.0:10) .* kJ)
energy_vec = SolarRequirement([collect(0.1:0.1:0.5); collect(1:5)] .* kJ)

# Set rates for birth and death
birth = fill(0.1/year, numSpecies)
death = fill(0.1/year, numSpecies)
longevity = 0.1
survival = 0.01
boost = 1.0
# Collect model parameters together
param = PopGrowth{typeof(unit(birth[1]))}(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(5.0km, 1e-9), numSpecies)
movement = BirthOnlyMovement(kernel, NoBoundary())

# Create species list, including their temperature preferences, seed abundance and native status
#opts = fill(274.0K, numSpecies)
opts = collect(270.0:279) .* K
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
p = 1 ./ energy_vec.energy 
abun = rand(Multinomial(individuals, p ./ sum(p)))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
abenv = tempgradAE(270.0K, 279.0K, grid, totalK, area, 0.0K/month)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create transition list
transitions = TransitionList(true)
addtransition!(transitions, UpdateEnergy(EcoSISTEM.update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for spp in eachindex(sppl.species.names)
    for loc in eachindex(abenv.habitat.matrix)
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, SeedDisperse(spp, loc))
    end
end
transitions = specialise_transition_list(transitions)

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)

# Simulation Parameters
burnin = 100years; times = 100years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@time simulate!(eco, burnin, timestep, specialise = true)
@time simulate_record!(abuns, eco, times, record_interval, timestep, specialise = true);

p = plot(sum(abuns[1, :, :], dims = 1)[1, :], grid = false, label = "", left_margin = 1.0*Plots.inch)
for i in 2:10
    p = plot!(sum(abuns[i, :, :], dims = 1)[1, :], grid = false, label = "", left_margin = 1.0*Plots.inch)
end
path = link_write!(handle, "PeatlandTest")
Plots.pdf(path)
finalise(handle)