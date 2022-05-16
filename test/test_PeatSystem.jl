using Peatland
using EcoSISTEM
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using Test

@testset "PeatSystem" begin
 
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
transitions = create_transition_list()
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
@test_nowarn eco = Ecosystem(sppl, abenv, rel, transitions = transitions, peatcache = true)
eco.cache.netmigration .+= 1
eco.cache.watermigration .+= 1.0m^3
Peatland.update_peat_environment!(eco, 1month)
@test sum(eco.cache.netmigration) == 0
@test sum(eco.cache.watermigration) == 0.0m^3
end