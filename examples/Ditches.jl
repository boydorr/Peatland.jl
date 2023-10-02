# using Revise
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
using Diversity

### HISTORIC RAINFALL DATA 2010 - 2020 (REPEATED) ###
function runPast(timestep::Unitful.Time, ditch = false; save = false, save_path = pwd())

    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

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
    height[end] = 5.0m # According to PlanAtt this is a more realistic height for Rhododendron
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
            addtransition!(transitions, WaterUse(spp, loc, 1.0, 0.02, 0.9/100.0mm))
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
        κ = 10 * 30.0m^2/month
        ν = 10 * 30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    for loc in neighbour_drains
        drainage = 0.0/month
        κ = 1 * 30.0m^2/month
        ν = 1 * 30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    for loc in not_drains
        if loc ∈ peat_locs
            κ = 10*30.0m^2/month
            ν = 10*30.0m^2/month
            fmax = 6.0/month
            kₛ = 0.09/100mm
            W0 = 0.6
            k2 = 5.0
        else
            κ = 10*30.0m^2/month
            ν = 10*30.0m^2/month
            fmax = 6.0/month
            kₛ = 0.09/100mm
            W0 = 0.5
            k2 = 5.0
        end
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        addtransition!(transitions, WaterFlux(loc, fmax, kₛ, W0, k2))
    end

    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 30year; times = 70year; 
    record_interval = 3months
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true);
    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    return abuns
end

## WITHOUT DITCHES ##
timestep = 1month
abuns = runPast(timestep);
@save "/home/claireh/sdc/Peatland/Peatland_past_noditch_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [7e6], ["A"])
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_new_test.pdf")

abuns_m = reshape(sum(abuns[mosses, :, end], dims = 1)[1, :], size(active))
heatmap(abuns_m', clim = (0,100))
abuns_s = reshape(sum(abuns[shrubs, :, end], dims = 1)[1, :], size(active))
heatmap(abuns_s', clim = (0,100))
abuns_t = reshape(sum(abuns[trees, :, end], dims = 1)[1, :], size(active))
heatmap(abuns_t')


## WITH DITCHES ##
abuns = runPast(timestep, true);
@save "/home/claireh/sdc/Peatland/Peatland_past_ditch_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [3.6e6], ["B"])
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_ditch_new.pdf")

### FUTURE RAINFALL SCENARIO 2010 - 2080 ###

function runFuture(timestep::Unitful.Time, ditch::Bool = false; save = false, save_path = pwd())

    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

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
    height[end] = 5.0m # According to PlanAtt this is a more realistic height for Rhododendron
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
            addtransition!(transitions, WaterUse(spp, loc, 1.0, 0.02, 0.9/100.0mm))
        end
    end

    # Add in location based transitions and ditches
    drains = [ditch_locs; river_locs ...]
    not_drains = setdiff(eachindex(abenv.habitat.h1.matrix), drains)
    for loc in drains
        drainage = 1.0/month
        κ = 10 * 30.0m^2/month
        ν = 10 * 30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν, ditch = ditch))
        addtransition!(transitions, Drainage(loc, drainage))
    end
    for loc in not_drains
        if loc ∈ peat_locs
            κ = 10*30.0m^2/month
            ν = 10*30.0m^2/month
            fmax = 6.0/month
            kₛ = 0.09/100mm
            W0 = 0.5
            k2 = 5.0
        else
            κ = 10*30.0m^2/month
            ν = 10*30.0m^2/month
            fmax = 6.0/month
            kₛ = 0.09/100mm
            W0 = 0.5
            k2 = 5.0
        end
        addtransition!(transitions, LateralFlow(loc, κ, ν, ditch = ditch))
        addtransition!(transitions, WaterFlux(loc, fmax, kₛ, W0, k2))
    end

    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)
    # Run simulation
    # Simulation Parameters
    burnin = 30year; times = 70year;
    record_interval = 3month
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true)

    @load "data/RainfallBudget_future.jld2"
    eco.abenv.budget = bud
    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);

    return abuns
end

## WITHOUT DITCHES ##
timestep = 1month
abuns = runFuture(timestep);
@save "/home/claireh/sdc/Peatland/Peatland_future_noditch_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [7e6], ["A"])
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_future_new.pdf")


## WITH DITCHES ##
abuns = runFuture(timestep, true);
@save "/home/claireh/sdc/Peatland/Peatland_future_ditch_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [3.6e6], ["B"])
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_future_ditch_new.pdf")


# Plotting for paper
# Heatmap for mosses
heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 300), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
@load "data/Peatland_past_noditch_new_test.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A", xlab = "Without drainage", ylab = "2010 - 2020 (repeated)", guide_position = :top)
@load "data/Peatland_past_ditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B", xlab = "With drainage", guide_position = :top)
@load "data/Peatland_future_noditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C", ylab = "2010 - 2080")
@load "data/Peatland_future_ditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Mosses_total_new.pdf")

# Heatmap for shrubs
heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 50), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
@load "data/Peatland_past_noditch_new_test.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A", xlab = "Without drainage", ylab = "2010 - 2020 (repeated)", guide_position = :top)
@load "data/Peatland_past_ditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B", xlab = "With drainage", guide_position = :top)
@load "data/Peatland_future_noditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C", ylab = "2010 - 2080")
@load "data/Peatland_future_ditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Others_total_new.pdf")

# Change over time
@load "data/Peatland_future_noditch_new.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
sumabuns2 = Float64.(reshape(sum(abuns[1:numMoss, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns2[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, (sumabuns .- sumabuns2)', c = :delta, clim = (-500, 500),
layout = (1, 2), size = (1200, 700), aspect_ratio = 1, grid = false, title = "A", titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)

sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
sumabuns2 = Float64.(reshape(sum(abuns[numMoss+1:end, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns2[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, (sumabuns .- sumabuns2)', c = :delta, clim = (-60, 60),
subplot = 2, aspect_ratio = 1, title = "B", titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
Plots.pdf("plots/Change_new.pdf")

