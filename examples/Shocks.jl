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
### FUTURE RAINFALL SCENARIO 2010 - 2080 plus summer drought in year 2060 ###

function runFuture(timestep::Unitful.Time, ditch::Bool = false; save = false, save_path = pwd())

    # Load data sets for angiosperms and mosses
    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)
    
    # Use outline raster of Cors Fochno to determine which grid cells are active
    file = "data/CF_outline.tif"
    cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    active = Array(.!isnan.(Array(cf)))
    #heatmap(active)

    # Set up parameters - grid and area aren't necessary but they are useful to know
    grid = size(cf); individuals=10_000_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid);
    numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub

    # Find all species in each category in order - for plotting
    mosses = 1:numMoss
    shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
    trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss

    # Extract size info from plant trait tables
    height = peat_spp[!, :Plant_height] .* m
    height[end] = 5.0m # According to PlanAtt this is a more realistic height for Rhododendron
    moss_height = moss_spp[!, :Len] .* mm

    #Set up how much energy each species consumes
    req1 = moss_height .* rand(TruncatedNormal(1.0, 0.1, 0, Inf), numMoss) .* mm ./m
    req2 = height .* rand(TruncatedNormal(10.0, 0.1, 0, Inf), numShrub) .* mm ./m
    energy_vec = WaterRequirement([req1; req2])

    # Set rates for birth and death and how much they impact survival and reproduction
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

    # Create species list, including their water level preferences, seed abundance and native status
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

    # Create abiotic environment from Topographic Wetness index, elevation and land cover category
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

    # If there are ditches, make sure they are lower than everything else around and start off with no water
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
    # Find neighbouring grid squares for all ditch locations
    xys = convert_coords.(drains, size(active, 1))
    neighbours = [Peatland.get_neighbours(active, i[1], i[2], 2) for i in xys]
    neighbours = unique(vcat(neighbours...))
    neighbour_drains = [convert_coords(n[1], n[2], size(active, 1)) for n in neighbours]
    not_drains = setdiff(eachindex(abenv.habitat.h1.matrix), vcat(drains, neighbour_drains))
    # Drainage ditches are drained and also act as sinks to the neighbour squares
    for loc in drains
        drainage = 1.0/month
        κ = 1 * 30.0m^2/month
        ν = 1 * 30.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    # Neighbour squares have their water pulled into the ditch
    for loc in neighbour_drains
        drainage = 0.0/month
        κ = 1 * 3.0m^2/month
        ν = 1 * 3.0m^2/month
        addtransition!(transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(transitions, Drainage(loc, drainage, κ, ν))
    end
    # Non-ditch squares act as normal
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
        # But also experience a certain amount of drying
        if ditch
            addtransition!(transitions, Dry(loc, 0.3, 1year))
        end
    end

    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 30year; times = 70year; timestep = 1month;
    record_interval = 3month
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true)

    @load "data/RainfallBudget_future.jld2"
    eco.abenv.budget = bud
    bud.matrix[:, :, 721:732] .*= 0.1

    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);

    return abuns
end

## WITHOUT DITCHES ##
timestep = 1month
abuns = runFuture(timestep);
@save "/home/claireh/sdc/Peatland/Peatland_future_shock_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [7e6], ["A"])
vline!([60], subplot = 1, color = :black, lty = :dot, label = "")
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
vline!([60], subplot = 2, color = :black, lty = :dot, label = "")
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
vline!([60], subplot = 3, color = :black, lty = :dot, label = "")
Plots.pdf("Abuns_future_shock_new.pdf")


### Past RAINFALL SCENARIO 2010 - 2080 plus summer drought in year 2060 ###

function runPast(ditch::Bool = false; save = false, save_path = pwd())
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
    # Add in species based transitions
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
            drainage = 0.001/month
        # end
        κ = 0.01/month
        λ = 0.01/month

        addtransition!(transitions, LateralFlow(loc, κ, λ))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 30year; times = 70year; timestep = 1month;
    record_interval = 3month
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    abuns1 = generate_storage(eco, lensim, 1)[:, :, :, 1]
    abuns2 = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true)

    fullbud = Array{typeof(bud.matrix[1]), 3}(undef, size(bud.matrix, 1), size(bud.matrix, 2), size(bud.matrix, 3) * 7)
    for i in 1:7
        fullbud[:, :, (i-1)*120+1:(i*120)] .= bud.matrix
    end
    bud.matrix = fullbud
    bud.matrix[:, :, 725:728] .*= 0.1

    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);

    return abuns
end

abuns = runPast();
@save "/home/claireh/sdc/Peatland/Peatland_past_shock_new.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [7e6], ["A"])
vline!([60], subplot = 1, color = :black, lty = :dot, label = "")
plot!(0:3/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
vline!([60], subplot = 2, color = :black, lty = :dot, label = "")
for i in trees
    plot!(0:3/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
vline!([60], subplot = 3, color = :black, lty = :dot, label = "")
Plots.pdf("Abuns_past_shock_new.pdf")

