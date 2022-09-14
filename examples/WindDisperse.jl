using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots

numSpecies = 1; grid = (10, 10); individuals=1000; area = 100.0*m^2; totalK = 10.0m^3/m^2

#Set up how much energy each species consumes
req1 = fill(0.1m^3, numSpecies)
energy_vec = VolWaterRequirement(req1)

#Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 100.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, NoBoundary())

# Create species list, including their temperature preferences, seed abundance and native status
pref1 = 100.0m^3
opts = fill(pref1, numSpecies)
vars = fill(10.0m^3, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
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
addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
for loc in eachindex(abenv.habitat.matrix)
    for spp in eachindex(sppl.species.names) 
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, WindDispersal(spp, loc, 1.0m, 0.3, 1.0m, 1.0m/s, 1.0m/s))
    end
end

lookup = PeatLookup(collect(1:100), collect(1:100), 100)

# Create ecosystem
eco = PeatSystem(sppl, abenv, rel, lookup, transitions = transitions)
eco.abundances.matrix[1, :] .= 0
eco.abundances.matrix[1, 1] = 100
# Run simulation
# Simulation Parameters
burnin = 1year; times = 5year; timestep = 1month;
record_interval = 1month; repeats = 1
lensim = length(0years:record_interval:times)
# Burnin
abuns = generate_storage(eco, lensim, 1)
@time simulate_record!(abuns, eco, times, record_interval, timestep)

plot(sum(abuns, dims = 2)[1, 1, :])
reshape_abuns = reshape(abuns, 1, 10, 10, lensim)
anim = @animate for i in 1:lensim
   heatmap(reshape_abuns[1, :, :, i], title = "Moss", clim = (0, maximum(abuns)))
end
gif(anim, "wind_disperse.gif", fps = 50)


using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Plots
using Profile
using ProfileView

function proftest(n::Int64)
    for i in 1:n
        numSpecies = 1; grid = (10, 10); individuals=1000; area = 100.0*m^2; totalK = 10.0m^3/m^2

        #Set up how much energy each species consumes
        req1 = fill(0.1m^3, numSpecies)
        energy_vec = VolWaterRequirement(req1)

        #Set rates for birth and death
        birth = 0.6/year
        death = 0.6/year
        longevity = 1.0
        survival = 0.1
        boost = 100.0
        # Collect model parameters together
        param = EqualPop(birth, death, longevity, survival, boost)

        # Create kernel for movement
        kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies)
        movement = BirthOnlyMovement(kernel, NoBoundary())

        # Create species list, including their temperature preferences, seed abundance and native status
        pref1 = 100.0m^3
        opts = fill(pref1, numSpecies)
        vars = fill(10.0m^3, numSpecies)
        traits = GaussTrait(opts, vars)
        native = fill(true, numSpecies)
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
        addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
        for loc in eachindex(abenv.habitat.matrix)
            for spp in eachindex(sppl.species.names) 
                addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
                addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
                addtransition!(transitions, GaussDispersal(spp, loc))
            end
        end

        # Create ecosystem
        eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
        eco.abundances.matrix[1, :] .= 0
        eco.abundances.matrix[1, 1] = 100
        # Run simulation
        # Simulation Parameters
        burnin = 1year; times = 5year; timestep = 1month;
        record_interval = 1month; repeats = 1
        lensim = length(0years:record_interval:times)
        abuns = generate_storage(eco, lensim, 1)
        @time simulate_record!(abuns, eco, times, record_interval, timestep);
    end
end

proftest(1)
ProfileView.@profview proftest(10)
