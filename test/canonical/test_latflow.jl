using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using AxisArrays
using Test

@testset "Diffusion" begin
    numSpecies = 2; grid = (10, 10); individuals=100;

    #Set up how much energy each species consumes
    req1 = 1.0mm; req2 = 10.0mm;
    energy_vec = WaterRequirement([req1, req2])

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
    pref1 = 10.0m^3; pref2 = 1.0m^3
    opts = [pref1, pref2]
    vars = [10.0m^3, 10.0m^3]
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    water_traits = GaussTrait(opts, vars)
    traits = TraitCollection2(water_traits, ele_traits)
    native = [true, false]
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    ele = AxisArray(fill(1.0, grid), Axis{:latitude}(1.0m:1.0m:10m), Axis{:longitude}(1.0m:1.0m:10m))
    abenv = peatlandAE(ele, 100.0m^3, 1000.0mm/m^2)
    abenv.habitat.h1.matrix[5, 5] = 90.0m^3

    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR2(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}())

    # Create new transition list
    transitions = TransitionList()
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
    drainage = 0.01/month
    κ = 0.001/month
    λ = 0.1/month
    for loc in eachindex(abenv.habitat.h1.matrix)
        for spp in eachindex(sppl.species.names) 
            addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(transitions, SeedDisperse(spp, loc))
            addtransition!(transitions, WaterUse(spp, loc, 0.001))
        end
        addtransition!(transitions, LateralFlow(loc, κ, λ))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    for i in 1:100
        display(heatmap(ustrip.(eco.abenv.habitat.h1.matrix), clim = (0, 100)))
        EcoSISTEM.update!(eco, 1month, transitions)
    end
    @test eco.abenv.habitat.h1.matrix[5,5] > 90.0m^3
end

@testset "Advection" begin
    numSpecies = 2; grid = (10, 10); individuals=100; 

    #Set up how much energy each species consumes
    req1 = 1.0mm; req2 = 10.0mm;
    energy_vec = WaterRequirement([req1, req2])

    #Set rates for birth and death
    birth = 0.1/year
    death = 0.1/year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies)
    movement = BirthOnlyMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = 10.0m^3; pref2 = 1.0m^3
    opts = [pref1, pref2]
    vars = [10.0m^3, 10.0m^3]
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    water_traits = GaussTrait(opts, vars)
    traits = TraitCollection2(water_traits, ele_traits)
    native = [true, false]
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    ele = AxisArray(reshape([i for i in 1.0:10 for j in 1:10], grid), Axis{:latitude}(1.0m:1.0m:10m), Axis{:longitude}(1.0m:1.0m:10m))
    abenv = peatlandAE(ele, 0.0m^3, 100.0mm/m^2)
    abenv.habitat.h1.matrix[:, 9] .= 200.0m^3
    
    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR2(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}())

    # Create new transition list
    transitions = TransitionList()
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))

    for loc in eachindex(abenv.habitat.h1.matrix)
        for spp in eachindex(sppl.species.names) 
            addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(transitions, SeedDisperse(spp, loc))
            addtransition!(transitions, WaterUse(spp, loc, 0.01))
        end
        drainage = 0.001/month
        κ = 0.01/month
        λ = 0.01/month

        addtransition!(transitions, LateralFlow(loc, κ, λ))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    for i in 1:100
        display(heatmap(ustrip.(eco.abenv.habitat.h1.matrix), clim = (0, 150)))
        EcoSISTEM.update!(eco, 1month, transitions)
    end

 
    @test eco.abenv.habitat.h1.matrix[:, 9] > eco.abenv.habitat.h1.matrix[:, 8] > eco.abenv.habitat.h1.matrix[:, 7] > eco.abenv.habitat.h1.matrix[:, 6]
end

@testset "Ditches" begin
    numSpecies = 2; grid = (10, 10); individuals=100; 

    #Set up how much energy each species consumes
    req1 = 1.0mm; req2 = 10.0mm;
    energy_vec = WaterRequirement([req1, req2])

    #Set rates for birth and death
    birth = 0.1/year
    death = 0.1/year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(10.0m, 10e-10), numSpecies)
    movement = BirthOnlyMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = 10.0m^3; pref2 = 1.0m^3
    opts = [pref1, pref2]
    vars = [10.0m^3, 10.0m^3]
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    water_traits = GaussTrait(opts, vars)
    traits = TraitCollection2(water_traits, ele_traits)
    native = [true, false]
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    ele = AxisArray(fill(2.0, grid), Axis{:latitude}(1.0m:1.0m:10m), Axis{:longitude}(1.0m:1.0m:10m))
    abenv = peatlandAE(ele, 100.0m^3, 100.0mm/m^2)
    abenv.habitat.h1.matrix[:, 5,] .= 10.0m^3
    abenv.habitat.h2.matrix[:, 5] .= 0.0

    ditches = findall(abenv.habitat.h1.matrix[1:end] .== 10.0m^3)

    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR2(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}())

    # Create new transition list
    transitions = TransitionList()
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))

    for loc in eachindex(abenv.habitat.h1.matrix)
        for spp in eachindex(sppl.species.names) 
            addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(transitions, SeedDisperse(spp, loc))
            addtransition!(transitions, WaterUse(spp, loc, 0.01))
        end
        if loc ∈ ditches
            drainage = 0.001/month
        else
            drainage = 0.001/month
        end
        κ = 0.01/month
        λ = 0.001/month

        addtransition!(transitions, LateralFlow(loc, κ, λ))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    for i in 1:100
        display(heatmap(ustrip.(eco.abenv.habitat.h1.matrix), clim = (0, 150)))
        EcoSISTEM.update!(eco, 1month, transitions)
    end

    @test eco.abenv.habitat.h1.matrix[:, 5] > eco.abenv.habitat.h1.matrix[:, 6]
    @test eco.abenv.habitat.h1.matrix[:, 6] < eco.abenv.habitat.h1.matrix[:, 7] < eco.abenv.habitat.h1.matrix[:, 8]
end