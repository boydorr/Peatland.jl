using Revise
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols
using Distributions
using DataFrames
using Plots
using JLD2
using Missings
using CSV
using Distances
using Shapefile
using Diversity

function buildEco(timestep::Unitful.Time; ditch = true)
    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

    moss_spp = filter(f -> (f.Len == findmax(moss_spp.Len)[1]) || (f.Len == findmin(moss_spp.Len)[1]), moss_spp)
    peat_spp = filter(f -> (f.Plant_height == findmax(peat_spp.Plant_height)[1]) || (f.Plant_height == findmin(peat_spp.Plant_height)[1]), peat_spp)

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
    req1 = moss_height .* rand(Normal(1, 0.1), numMoss) .* mm ./m
    req2 = height .* rand(Normal(10.0, 0.1), numShrub) .* mm ./m
    energy_vec = WaterRequirement([req1; req2])

    # Set rates for birth and death
    birth = 0.1/year
    death = 0.1/year
    longevity = 0.1
    survival = 0.01
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(20.0m, 1e-5), numSpecies)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = rand(TruncatedNormal(0.4, 0.01, 0.0, 1.0), numMoss)
    pref2 = rand(TruncatedNormal(0.3, 0.01, 0.0, 1.0), numShrub)
    opts = [pref1; pref2]
    vars = [fill(0.01, numMoss); fill(0.01, numShrub)]
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

    # Create abiotic environment - even grid of one temperature
    file = "data/CF_elevation_smooth.tif"
    ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    ele.data[ele.data .< 0] .= 0

    file = "data/CF_TPI_smooth.tif"
    tpi = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    tpi.data[tpi.data .< 0] .= 0
    tpi.data ./= maximum(tpi.data)

    @load "data/RainfallBudget_burnin.jld2"
    bud.matrix .*= (timestep / month)
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(tpi, ele, soil, bud, active)

    # If there are ditches, make sure they are lower than everything else around and filled with water
    if ditch
        file = "data/Ditches_full.tif"
        ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        ditch_locations = findall(.!isnan.(ditches))
        ditch_locs = [convert_coords(d[1], d[2], size(cf, 1)) for d in ditch_locations]
        file = "data/MainRivers.tif"
        rivers = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        river_locations = findall(.!isnan.(rivers))
        river_locs = [convert_coords(r[1], r[2], size(cf, 1)) for r in river_locations]
        locs = [ditch_locs; river_locs ...]
        # locs = ditch_locs
        abenv.habitat.h1.matrix[locs] .= 0.0
        abenv.habitat.h2.matrix[locs] .= 0
        abenv.active[ditch_locations] .= false
        abenv.active[river_locations] .= false
    else
        ditch_locs = []
        river_locs = []
    end
    # heatmap(ustrip.(abenv.habitat.h1.matrix)')

    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR3(Gauss{Float64}(), NoRelContinuous{Float64}(), soilmatch{Int64}())

    # Create new transition list
    transitions = TransitionList(true)
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
    active_squares = findall(active[1:end])
    peat_squares = findall(soil .== 11)
    peat_locs = [convert_coords(d[1], d[2], size(cf, 1)) for d in peat_squares]
    # Add in species based transitions (only for active squares)
    for spp in eachindex(sppl.species.names) 
        for loc in active_squares
            addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(transitions, SeedDisperse(spp, loc))
            if spp > numMoss
                addtransition!(transitions, Invasive(spp, loc, 10.0/28days))
            end
        end
    end
    # Water needs to be used everywhere (with a background rate for where we aren't modelling plants)
    for spp in eachindex(sppl.species.names) 
        for loc in eachindex(abenv.habitat.h1.matrix)
            addtransition!(transitions, WaterUse(spp, loc, 0.02, 0.001))
        end
    end

    # # Add in location based transitions and ditches
    drains = [ditch_locs; river_locs ...]
    not_drains = setdiff(eachindex(abenv.habitat.h1.matrix), drains)
    for loc in drains
        drainage = 1.0/month
        κ = 100*30.0m^2/month
        ν = 100*30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        addtransition!(transitions, Drainage(loc, drainage))
    end
    for loc in not_drains
        if loc ∈ peat_locs
            κ = 100.0m^2/month
            ν = 100.0m^2/month
            fmax = 1.0/month
            kₛ = 3e-4/m^3
            ϵ = 1.0m^3
        else
            κ = 100.0m^2/month
            ν = 100.0m^2/month
            fmax = 1.0/month
            kₛ = 2e-4/m^3
            ϵ = 1.0m^3
        end
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        addtransition!(transitions, WaterFlux(loc, fmax, kₛ, ϵ))
    end


    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    return eco
end
timestep = 1month
eco = buildEco(timestep, ditch = false)
envs = zeros(100)
for i in 1:100
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 10), subplot = 2)
plot(envs, ylim= (0,1))

abuns = reshape(sum(eco.abundances.matrix, dims = 1), size(active))
heatmap(abuns')

abuns_m = reshape(sum(eco.abundances.matrix[mosses, :], dims = 1)[1, :], size(active))
heatmap(abuns_m')
abuns_s = reshape(sum(eco.abundances.matrix[shrubs, :], dims = 1)[1, :], size(active))
heatmap(abuns_s')

eco = buildEco(timestep, ditch = true)
envs = zeros(20)
for i in 1:20
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 100), subplot = 2)
plot(envs)