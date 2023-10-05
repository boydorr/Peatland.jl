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
using Diversity

function buildEco(timestep::Unitful.Time; ditch = true, specialise=true)
    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

    # moss_spp = filter(f -> (f.Len == findmax(moss_spp.Len)[1]) || (f.Len == findmin(moss_spp.Len)[1]), moss_spp)
    # peat_spp = filter(f -> (f.Plant_height == findmax(peat_spp.Plant_height)[1]) || (f.Plant_height == findmin(peat_spp.Plant_height)[1]), peat_spp)
    species = rand(1:nrow(peat_spp), 10)
    peat_spp = peat_spp[species, :]
    species = rand(1:nrow(moss_spp), 10)
    moss_spp = moss_spp[species, :]
    file = "data/CF_outline.tif"
    cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    active = Array(.!isnan.(Array(cf)))
    #heatmap(active)

    grid = size(cf); individuals=10_000_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid);
    numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub

    mosses = 1:numMoss
    shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
    trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss

    height = peat_spp[!, :Plant_height] .* m
    moss_height = moss_spp[!, :Len] .* mm

    #Set up how much energy each species consumes
    req1 = moss_height .* rand(TruncatedNormal(1.0, 0.1, 0, Inf), numMoss) .* mm ./m
    req2 = height .* rand(TruncatedNormal(10.0, 0.1, 0, Inf), numShrub) .* mm ./m
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
    pref1 = rand(TruncatedNormal(0.6, 0.05, 0.0, 1.0), numMoss)
    pref2 = rand(TruncatedNormal(0.3, 0.05, 0.0, 1.0), numShrub)
    opts = [pref1; pref2]
    vars = [fill(0.05, numMoss); fill(0.05, numShrub)]
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
    bud.matrix .= mean(bud.matrix)
    # # Need to detrend and centre around mean
    # newbud = zeros(typeof(1.0mm), size(bud.matrix))
    # newbud[:, :, 1] .= bud.matrix[:, :, 1]
    # meanbud = mean(bud.matrix)
    # for i in 2:size(bud.matrix, 3)
    #     newbud[:, :, i] .= (bud.matrix[:, :, i] .- bud.matrix[:, :, i-1])
    #     newbud[:, :, i] .+= meanbud
    # end
    # bud.matrix .= abs.(newbud)
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(tpi, ele, soil, bud, active)

    # If there are ditches, make sure they are lower than everything else around and filled with water
    if ditch
        file = "data/Ditches_full.tif"
        ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        ditch_locations = findall(.!isnan.(ditches))
        ditch_locs = [convert_coords(d[1], d[2], size(active, 1)) for d in ditch_locations]
        file = "data/MainRivers.tif"
        rivers = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        river_locations = findall(.!isnan.(rivers))
        river_locs = [convert_coords(r[1], r[2], size(active, 1)) for r in river_locations]
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
            # if spp > numMoss
            addtransition!(transitions, Invasive(spp, loc, 1.0/30days))
            # end
        end
    end
    # Water needs to be used everywhere (with a background rate for where we aren't modelling plants)
    for spp in eachindex(sppl.species.names) 
        for loc in eachindex(abenv.habitat.h1.matrix)
            addtransition!(transitions, WaterUse(spp, loc, 1.0, 0.02, 0.6/100.0mm))
        end
    end

    # # Add in location based transitions and ditches
    drains = [ditch_locs; river_locs ...]
    xys = convert_coords.(drains, size(active, 1))
    neighbours = [Peatland.get_neighbours(active, i[1], i[2], 2) for i in xys]
    neighbours = unique(vcat(neighbours...))
    neighbour_drains = [convert_coords(n[1], n[2], size(active, 1)) for n in neighbours]
    not_drains = setdiff(eachindex(abenv.habitat.h1.matrix), vcat(drains, neighbour_drains))
    for loc in drains
        drainage = 1.0/month
        κ = 1 * 30.0m^2/month
        ν = 1 * 30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    for loc in neighbour_drains
        drainage = 0.0/month
        κ = 1 * 3.0m^2/month
        ν = 1 * 3.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    for loc in not_drains
        κ = 10*30.0m^2/month
        ν = 10*30.0m^2/month
        fmax = 6.0/month
        kₛ = 0.06/100mm
        W0 = 0.5
        k2 = 5.0
       
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        addtransition!(transitions, WaterFlux(loc, fmax, kₛ, W0, k2))
    end

    if specialise
        transitions = specialise_transition_list(transitions)
    end
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    return eco
end
timestep = 1month
eco = buildEco(timestep, ditch = false)
envs = zeros(120)
for i in 1:120
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
envs = zeros(120)
for i in 1:120
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 100), subplot = 2)
plot(envs)


eco = buildEco(timestep, specialise = false)
envs = zeros(36)
for i in 1:12
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 100), subplot = 2)

for loc in active_squares
    addtransition!(eco.transitions, Peatland.Dry(loc, 0.05, 1year))
end
for i in 13:24
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 100), subplot = 2)

for loc in active_squares
    addtransition!(eco.transitions, Peatland.Rewet(loc, 0.05, 1year, 1.0))
end
for i in 25:36
    EcoSISTEM.update!(eco, timestep, eco.transitions)
    envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
end
heatmap(ustrip.(eco.abenv.habitat.h1.matrix)', clim = (0, 1), size = (1000, 600), layout = 2)
heatmap!(ustrip.(eco.cache.surfacewater)', clim = (0, 100), subplot = 2)
plot(envs)