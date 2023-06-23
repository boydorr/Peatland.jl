using Revise
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
using Distances
using Shapefile
using Diversity

### HISTORIC RAINFALL DATA 2010 - 2020 (REPEATED) ###
function runPast(paramDict::Dict, ditch = false; reps = 10, save = false, save_path = pwd())
    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

    file = "data/CF_outline.tif"
    cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    active = Array(.!isnan.(Array(cf)))
    #heatmap(active)

    grid = size(cf); individuals=1_000_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid);
    numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub

    mosses = 1:numMoss
    shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
    trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss

    height = peat_spp[!, :Plant_height] .* m
    moss_height = moss_spp[!, :Len] .* mm

    #Set up how much energy each species consumes
    req1 = moss_height .* rand(Normal(1.0, 0.1), numMoss) .* mm ./m
    req2 = height .* rand(Normal(10.0, 0.1), numShrub) .* mm ./m
    energy_vec = WaterRequirement([req1; req2])


    # Set rates for birth and death
    birth = paramDict.birth
    death = paramDict.death
    longevity = paramDict.longevity
    survival = paramDict.survival
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(paramDict.move, 1e-5), numSpecies)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = rand(Normal(90.0, 13.0), numMoss) .* m^3 
    pref2 = rand(Normal(70.0, 18.0), numShrub) .* m^3 
    opts = [pref1; pref2]
    vars = [fill(10.0m^3, numMoss); fill(10.0m^3, numShrub)]
    water_traits = GaussTrait(opts, vars)
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    soilDict = Dict("hygrophilous" => [8, 11], "terrestrial" => [1, 4, 5], "terrestrial/hygrophilous" => [1, 4, 5, 8, 11])
    soil_pref1 = fill([8, 11], numMoss); soil_pref2 = [soilDict[x] for x in peat_spp[!, :Habitat]]
    soil_pref = soiltrait([soil_pref1; soil_pref2])
    traits = TraitCollection3(water_traits, ele_traits, soil_pref)
    native = fill(true, numSpecies)
    p = 1 ./ energy_vec.energy 
    abun = rand(Multinomial(individuals, p ./ sum(p)))
    # abun = [fill(div(individuals, numMoss), numMoss); fill(0, numShrub)]
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    for i in reps
        # Create abiotic environment - even grid of one temperature
        file = "data/CF_elevation_smooth.tif"
        ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
        ele.data[ele.data .< 0] .= 0

        file = "data/CF_TPI_smooth.tif"
        tpi = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
        tpi.data[tpi.data .< 0] .= 0

        @load "data/RainfallBudget_burnin.jld2"
        file = "data/LCM.tif"
        soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
        soil = Int.(soil)
        abenv = peatlandAE(ele, soil, 10.0m^3, bud, active) 
        abenv.habitat.h1.matrix .*= tpi.data
        abenv.habitat.h1.matrix[abenv.habitat.h1.matrix .< 10.0m^3] .= 10.0m^3

        # If there are ditches, make sure they are lower than everything else around and filled with water
        if ditch
            file = "data/Ditches_full.tif"
            ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
            ditch_locations = findall(.!isnan.(ditches))
            locs = [convert_coords(d[1], d[2], size(cf, 1)) for d in ditch_locations]
            abenv.habitat.h1.matrix[locs] .= 0.0m^3
            abenv.habitat.h2.matrix[locs] .= 0
            abenv.active[ditch_locations] .= false
        else
            ditch_locations = []
        end
        # heatmap(ustrip.(abenv.habitat.h1.matrix)')

        # Set relationship between species and environment (gaussian)
        rel = multiplicativeTR3(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}(), soilmatch{Int64}())

        # Create new transition list
        transitions = TransitionList(true)
        addtransition!(transitions, UpdateEnergy(update_energy_usage!))
        addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
        active_squares = findall(active[1:end])
        peat_squares = findall(soil .== 11)
        # Add in species based transitions (only for active squares)
        for spp in eachindex(sppl.species.names) 
            for loc in active_squares
                addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
                addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
                addtransition!(transitions, SeedDisperse(spp, loc))
                addtransition!(transitions, WaterUse(spp, loc, 0.01))
                if spp > numMoss
                    addtransition!(transitions, Invasive(spp, loc, 10.0/28days))
                end
            end
        end
        # Add in location based transitions and ditches
        for loc in eachindex(abenv.habitat.h1.matrix)
            # if loc ∈ ditch_locations
            #     drainage = 1.0/month
            # else
                drainage = 0.029/month
            # end
            κ = 0.01m^2/month
            λ = 0.01m^2/month

            addtransition!(transitions, LateralFlow(loc, κ, λ))
            addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
        end

        transitions = specialise_transition_list(transitions)
        # Create ecosystem
        eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

        # Run simulation
        # Simulation Parameters
        burnin = 30year; times = 70year; timestep = 1year;
        record_interval = 1year
        lensim = length(0years:record_interval:times)
        # Burnin
        abuns = generate_storage(eco, lensim, 10)[:, :, :, 1]
        @time simulate!(eco, burnin, timestep, specialise = true);
        @time simulate_record!(abuns[:, :, i], eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    end
    return abuns
end

## WITHOUT DITCHES ##
for i in 1:10
    paramDict = Dict()
    abuns = runPast();