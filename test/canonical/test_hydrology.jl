using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using Test

@testset "Hydrology" begin

    grid = (75,64); individuals=5_000; area = 12.0km^2; totalK = 0.01mm/m^2
    numMoss = 5; numShrub = 28; numSpecies = numMoss + numShrub

    shrub_size = fill(5.0m^2, numShrub)
    moss_height = fill(1.0, numMoss)
    moss_size = moss_height .* 1.0m^2
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
    abenv = peatlandAE(100.0m^3, grid, totalK, area)

    # Set relationship between species and environment (gaussian)
    rel = Gauss{typeof(1.0m^3)}()

    # Create new transition list
    transitions = TransitionList()
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
    active_squares = 1:prod(grid)
    for loc in active_squares
        for spp in eachindex(sppl.species.names) 
            addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(transitions, SeedDisperse(spp, loc))
            addtransition!(transitions, WaterUse(spp, loc, 0.01))
            if spp > numMoss
                addtransition!(transitions, Invasive(spp, loc, 10.0/28days))
            end
        end
        addtransition!(transitions, LateralFlow(abenv, loc, 1/600days))
        addtransition!(transitions, WaterFlux(loc, 0.625/month))
    end


    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 10year; times1 = 5year; times2 = 5year; timestep = 1months;
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

    sum_hydro = sum(eco.abenv.habitat.matrix)
    for loc in active_squares
        addtransition!(transitions, Peatland.Dry(loc, 0.5, 1year))
    end

    @time simulate_record!(abuns1, eco, times1, record_interval, timestep);
    sum_dried = sum(eco.abenv.habitat.matrix)
    @test sum_dried < sum_hydro

    for loc in active_squares
        addtransition!(transitions, Peatland.Rewet(loc, 0.5, 1year, 150.0m^3))
    end

    @time simulate_record!(abuns2, eco, times2, record_interval, timestep);
    @test sum(eco.abenv.habitat.matrix) > sum_dried

    # Test for Moss vs Shrub numbers
    @test sum(abuns[1:numMoss, :, end]) > sum(abuns1[1:numMoss, :, end])
    @test sum(abuns1[1:numMoss, :, end]) < sum(abuns2[1:numMoss, :, end])

    @test sum(abuns[(numMoss+1):end, :, end]) < sum(abuns1[(numMoss+1):end, :, end])
    @test sum(abuns1[(numMoss+1):end, :, end]) > sum(abuns2[(numMoss+1):end, :, end])
end