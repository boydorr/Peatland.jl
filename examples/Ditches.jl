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
function runPast(ditch = false; save = false, save_path = pwd())
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

    height = peat_spp[!, :Plant_height]
    shrub_size = height .* 1.0m^2
    moss_height = moss_spp[!, :Len] .* mm
    moss_size = uconvert.(m^2, moss_height .* 0.1m)

    #Set up how much energy each species consumes
    req1 = moss_size .* rand(Normal(10.0, 0.1), numMoss) .* mm ./m^2
    req2 = shrub_size .* rand(Normal(10.0, 0.1), numShrub) .* mm ./m^2
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

    envs = zeros(10)
    for i in 1:10
        EcoSISTEM.update!(eco, timestep, transitions)
        envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
    end
    heatmap(ustrip.(eco.abenv.habitat.h1.matrix)')
    #plot(envs)

    # hab = ustrip.(eco.abenv.habitat.h1.matrix)
    # #hab[.!active] .= NaN
    # heatmap(hab', clim = (0, 100))
    # Plots.pdf("plots/Hab.pdf")

    # hab = ustrip.(eco.abenv.habitat.h2.matrix)
    # hab[.!active] .= NaN
    # heatmap(hab)
    # Plots.pdf("plots/Ele.pdf")

    # Run simulation
    # Simulation Parameters
    burnin = 20year; times = 70year; timestep = 1month;
    record_interval = 3months
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true);
    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    return abuns
end

## WITHOUT DITCHES ##
abuns = runPast();
@save "/home/claireh/sdc/Peatland/Peatland_past_noditch.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns.pdf")

## WITH DITCHES ##
abuns = runPast(true);
@save "/home/claireh/sdc/Peatland/Peatland_past_ditch.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_ditch.pdf")

### FUTURE RAINFALL SCENARIO 2010 - 2080 ###

function runFuture(ditch::Bool = false; save = false, save_path = pwd())
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

    height = peat_spp[!, :Plant_height]
    shrub_size = height .* 1.0m^2
    moss_height = moss_spp[!, :Len] .* mm
    moss_size = uconvert.(m^2, moss_height .* 0.1m)

    #Set up how much energy each species consumes
    req1 = moss_size .* rand(Normal(10.0, 0.1), numMoss) .* mm ./m^2
    req2 = shrub_size .* rand(Normal(10.0, 0.1), numShrub) .* mm ./m^2
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
    file = "data/CF_elevation_full.tif"
    ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    ele.data[ele.data .< 0] .= 0
    # If there are ditches, make sure they are lower than everything else around
    if ditch
        file = "data/Ditches_full.tif"
        ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        ditch_locations = findall(.!isnan.(ditches))
        locs = [convert_coords(d[1], d[2], size(cf, 1)) for d in ditch_locations]
        ele[locs] .= 0
    end
    file = "data/CF_TPI_corrected.tif"
    tpi = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    tpi.data[tpi.data .< 0] .= 0

    @load "data/RainfallBudget_burnin.jld2"
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(ele, soil, 10.0m^3, bud, active) 
    abenv.habitat.h1.matrix .*= tpi.data
    abenv.habitat.h1.matrix[abenv.habitat.h1.matrix .< 10.0m^3] .= 10.0m^3
    # heatmap(ustrip.(abenv.habitat.h1.matrix)')

    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR3(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}(), soilmatch{Int64}())

    # Create new transition list
    transitions = TransitionList(true)
    addtransition!(transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(transitions, UpdateEnvironment(update_peat_environment!))
    active_squares = findall(active[1:end])
    peat_squares = findall(soil .== 1)
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
    # Add in location based transitions
    for loc in eachindex(abenv.habitat.h1.matrix)
        if loc ∈ peat_squares
            drainage = 0.00001/month; flow = 0.05/month
        else
            drainage = 0.006/month; flow = 0.5/month
        end
        addtransition!(transitions, LateralFlow(abenv, abenv.habitat.h2, loc, flow))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    # Add ditch drainage
    if ditch
        file = "data/Ditches_full.tif"
        ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
        ditch_locations = findall(.!isnan.(ditches))
        # Alter drainage in these cells
        for d in ditch_locations
            loc = convert_coords(d[1], d[2], size(cf, 1))
            flow = transitions.setup[loc*2]
            flux = transitions.setup[loc*2+1]
            flow.prob = 1.0/month
            flux.prob = 1.0/month
        end
    end
    transitions = specialise_transition_list(transitions)
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 20year; times = 70year; timestep = 1month;
    record_interval = 1month
    lensim = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    abuns1 = generate_storage(eco, lensim, 1)[:, :, :, 1]
    abuns2 = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true)

    @load "data/RainfallBudget_future.jld2"
    eco.abenv.budget = bud
    @time simulate_record!(abuns, eco, times, record_interval, timestep, save = save, save_path = save_path, specialise = true);

    return abuns
end

## WITHOUT DITCHES ##
abuns = runFuture();
@save "/home/claireh/sdc/Peatland/Peatland_future_noditch.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_future.pdf")


## WITH DITCHES ##
abuns = runFuture(true);
@save "/home/claireh/sdc/Peatland/Peatland_future_ditch.jld2" abuns=abuns[:, :, [12,end]]

plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[mosses, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
for i in shrubs
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[shrubs, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
for i in trees
    plot!(0:1/12:70,sum(abuns[i, :, :], dims = 1)[1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:1/12:70, mean(sum(abuns[trees, :, :], dims = 2)[:, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
Plots.pdf("Abuns_future_ditch.pdf")


# Plotting for paper
# Heatmap for mosses
heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 600), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
@load "data/Peatland_past_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A", xlab = "Without drainage", ylab = "2010 - 2020 (repeated)", guide_position = :top)
@load "data/Peatland_past_ditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B", xlab = "With drainage", guide_position = :top)
@load "data/Peatland_future_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C", ylab = "2010 - 2080")
@load "data/Peatland_future_ditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Mosses_total.pdf")

# Heatmap for shrubs
heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 50), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
@load "data/Peatland_past_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A", xlab = "Without drainage", ylab = "2010 - 2020 (repeated)", guide_position = :top)
@load "data/Peatland_past_ditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B", xlab = "With drainage", guide_position = :top)
@load "data/Peatland_future_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C", ylab = "2010 - 2080")
@load "data/Peatland_future_ditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Others_total.pdf")

# Change over time
@load "data/Peatland_future_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
sumabuns2 = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns2[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, (sumabuns .- sumabuns2)', c = :delta, clim = (-500, 500),
layout = (1, 2), size = (1200, 700), aspect_ratio = 1, grid = false)

sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, 1], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
sumabuns2 = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns2[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, (sumabuns .- sumabuns2)', c = :delta, clim = (-60, 60),
subplot = 2, aspect_ratio = 1)
Plots.pdf("plots/Change.pdf")


@load "data/Peatland_future_noditch.jld2"
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
@load "data/Peatland_future_ditch.jld2"
sumabuns2 = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns2[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, (sumabuns .- sumabuns2)', c = :delta, clim = (-200, 200),
layout = (1, 2), size = (1200, 700), aspect_ratio = 1, grid = false)