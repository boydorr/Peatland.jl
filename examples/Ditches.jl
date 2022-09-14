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
function runsim()
    JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
    JLD2.@load("data/Peat_30_moss.jld2", moss_spp)

    file = "data/CF_outline.tif"
    cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    active = Array(.!isnan.(Array(cf)))
    #heatmap(active)

    grid = size(cf); individuals=100_000; area = step(cf.axes[1]) * step(cf.axes[2]) * prod(grid);
    numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub

    mosses = 1:numMoss
    shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
    trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss

    height = peat_spp[!, :Plant_height]
    shrub_size = height .* 1.0m^2
    moss_height = moss_spp[!, :Len] .* mm
    moss_size = uconvert.(m^2, moss_height .* 0.1m)

    #Set up how much energy each species consumes
    req1 = moss_size .* rand(Normal(100.0, 0.1), numMoss) .* mm ./m^2
    req2 = shrub_size .* rand(Normal(1.0, 0.1), numShrub) .* mm ./m^2
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
    kernel = fill(GaussianKernel(10.0m, 1e-3), numSpecies)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = rand(Normal(120.0, 10.0), numMoss) .* m^3 
    pref2 = abs.(rand(Normal(60.0, 20.0), numShrub)) .* m^3 
    opts = [pref1; pref2]
    vars = [fill(5.0m^3, numMoss); fill(10.0m^3, numShrub)]
    water_traits = GaussTrait(opts, vars)
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    soilDict = Dict("hygrophilous" => [8, 11], "terrestrial" => [1, 4, 5], "terrestrial/hygrophilous" => [1, 4, 5, 8, 11])
    soil_pref1 = fill([8, 11], numMoss); soil_pref2 = [soilDict[x] for x in peat_spp[!, :Habitat]]
    soil_pref = soiltrait([soil_pref1; soil_pref2])
    traits = TraitCollection3(water_traits, ele_traits, soil_pref)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    # abun = [fill(div(individuals, numMoss), numMoss); fill(0, numShrub)]
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    file = "data/CF_elevation_full.tif"
    ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    ele.data[ele.data .< 0] .= 0

    file = "data/CF_TPI_corrected.tif"
    tpi = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    tpi.data[tpi.data .< 0] .= 0

    @load "data/RainfallBudget.jld2"
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(ele, soil, 10.0m^3, bud, active) 
    abenv.habitat.h1.matrix .*= tpi.data
    abenv.habitat.h1.matrix[abenv.habitat.h1.matrix .< 10.0m^3] .= 10.0m^3
    heatmap(ustrip.(abenv.habitat.h1.matrix)')

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
            drainage = 0.00001/month; flow = 0.0001/month
        else
            drainage = 0.006/month; flow = 0.5/month
        end
        addtransition!(transitions, LateralFlow(abenv, abenv.habitat.h2, loc, flow))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end
    transitions = specialise_transition_list(transitions)
    # Add ditch drainage
    # file = "data/Ditches.tif"
    # ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
    # ditch_locations = findall(.!isnan.(ditches))
    # for d in ditch_locations
    #     loc = convert_coords(d[2], d[1], size(cf, 1))
    #     addtransition!(transitions, Dry(loc, 1.0, 1month))
    # end
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # envs = zeros(lensim)
    # for i in 1:lensim
    #     EcoSISTEM.update!(eco, timestep, transitions)
    #     envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
    # end
    # plot(envs)

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
    burnin = 5year; times = 10year; timestep = 1month;
    record_interval = 1month
    lensim = length(0years:record_interval:burnin)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep, specialise = true);
    @time simulate_record!(abuns, eco, burnin, record_interval, timestep, save = false);
    return abuns
end

abuns = runsim()

# JLD2.save("data/CF_abuns.jld2", "abuns", abuns)
# JLD2.@load "data/CF_abuns.jld2"

sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', size = (1000, 800), aspect_ratio = 1,
margin = 1.0*Plots.inch)
Plots.pdf("Mosses.pdf")

sumabuns = Float64.(reshape(sum(abuns[numMoss+1:end, :, end], dims = 1), size(abenv.habitat.h1.matrix)))
sumabuns[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', size = (1000, 800), aspect_ratio = 1,
margin = 1.0*Plots.inch)
Plots.pdf("Others.pdf")

plot(0:1/12:5,sum(abuns[1:numMoss, :, :, :, 1], dims = (1,2))[1, 1, :], grid = false, label = "",
layout = (2, 1), title = "Moss", left_margin = 1.0*Plots.inch)
plot!(0:1/12:5, sum(abuns[numMoss+1:end, :, :, :, 1], dims = (1,2))[1, 1, :], label = "",
subplot = 2, color = 2, grid = false, title = "Other")
Plots.pdf("Abuns.pdf")


# cf_outline = Shapefile.shapes(Shapefile.Table("data/CorsFochno.shp"))
# ys = ustrip.(cf.axes[1].val); xs = ustrip.(cf.axes[2].val)
# heatmap(ys, xs, Shrub[:, :, end]', layout = (@layout [a b]), size = (1000, 400))
# plot!(cf_outline, fillcolor = false, linecolor = :white, subplot = 1)
# heatmap!(ys, xs, Moss[:, :, end]', subplot =2)
# plot!(cf_outline, fillcolor = false, linecolor = :white, subplot =2)


function runsim()
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
    birth = 0.6/year
    death = 0.6/year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(10.0m, 1e-3), numSpecies)
    movement = BirthOnlyMovement(kernel, NoBoundary())

    # Create species list, including their temperature preferences, seed abundance and native status
    pref1 = rand(Normal(120.0, 0.1), numMoss) .* m^3 
    pref2 = abs.(rand(Normal(70.0, 20.0), numShrub)) .* m^3 
    opts = [pref1; pref2]
    vars = [fill(5.0m^3, numMoss); fill(10.0m^3, numShrub)]
    water_traits = GaussTrait(opts, vars)
    ele_traits = GaussTrait(fill(1.0, numSpecies), fill(20.0, numSpecies))
    soilDict = Dict("hygrophilous" => [8, 11], "terrestrial" => [1, 4, 5], "terrestrial/hygrophilous" => [1, 4, 5, 8, 11])
    soil_pref1 = fill([8, 11], numMoss); soil_pref2 = [soilDict[x] for x in peat_spp[!, :Habitat]]
    soil_pref = soiltrait([soil_pref1; soil_pref2])
    traits = TraitCollection3(water_traits, ele_traits, soil_pref)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    # abun = [fill(div(individuals, numMoss), numMoss); fill(0, numShrub)]
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    file = "data/CF_elevation_full.tif"
    ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    ele.data[ele.data .< 0] .= 0

    file = "data/CF_TPI_corrected.tif"
    tpi = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    tpi.data[tpi.data .< 0] .= 0

    @load "data/RainfallBudget.jld2"
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(ele, soil, 10.0m^3, bud, active) 
    abenv.habitat.h1.matrix .*= tpi.data
    abenv.habitat.h1.matrix[abenv.habitat.h1.matrix .< 10.0m^3] .= 10.0m^3
    heatmap(ustrip.(abenv.habitat.h1.matrix)')

    # Set relationship between species and environment (gaussian)
    rel = multiplicativeTR3(Gauss{typeof(1.0m^3)}(), NoRelContinuous{Float64}(), soilmatch{Int64}())

    # Create new transition list
    transitions = TransitionList()
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
            drainage = 0.00001/month; flow = 0.0001/month
        else
            drainage = 0.006/month; flow = 0.5/month
        end
        addtransition!(transitions, LateralFlow(abenv, abenv.habitat.h2, loc, flow))
        addtransition!(transitions, WaterFlux(loc, drainage, 150.0m^3))
    end

    # Add ditch drainage
    # file = "data/Ditches.tif"
    # ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
    # ditch_locations = findall(.!isnan.(ditches))
    # for d in ditch_locations
    #     loc = convert_coords(d[2], d[1], size(cf, 1))
    #     addtransition!(transitions, Dry(loc, 1.0, 1month))
    # end
    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # envs = zeros(lensim)
    # for i in 1:lensim
    #     EcoSISTEM.update!(eco, timestep, transitions)
    #     envs[i] = ustrip.(mean(eco.abenv.habitat.h1.matrix[eco.abenv.active]))
    # end
    # plot(envs)

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
    burnin = 5year; times = 5year; timestep = 1month;
    record_interval = 1month
    lensim = length(0years:record_interval:burnin)
    lensim2 = length(0years:record_interval:times)
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    abuns1 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
    abuns2 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
    @time simulate!(eco, burnin, timestep)
    @time simulate_record!(abuns, eco, burnin, record_interval, timestep, save = false);
    for loc in active_squares
        addtransition!(transitions, Peatland.Dry(loc, 0.01, 1year))
    end
    @time simulate_record!(abuns1, eco, times, record_interval, timestep, save = false);
    for loc in active_squares
        addtransition!(transitions, Peatland.Rewet(loc, 0.01, 1year, 150.0m^3))
    end
    @time simulate_record!(abuns2, eco, times, record_interval, timestep, save = false);
    abuns = cat(abuns, abuns1, abuns2,  dims = 3)
    abuns = reshape(abuns, (numSpecies, grid[1], grid[2], 61*3))
    
    return abuns
end

abuns = runsim()



times = 1:1/12:(15+2/12)
p = plot(times, sum(abuns[1:numMoss, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :], grid = false, label = "",
layout = (3, 1), title = "Moss", size = (1000, 1200), margin = 1.0*Plots.inch)
p = plot!(times, sum(abuns[shrubs, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :], label = "",
subplot = 2, color = 2, grid = false, title = "Shrubs, Herbs and Grasses", margin = 1.0*Plots.inch)
p = plot!(times, sum(abuns[trees, :, :, 13:end, 1], dims = (1,2,3))[1, 1, 1, :], label = "",
subplot = 3, color = 3, grid = false, title = "Trees", margin = 1.0*Plots.inch)
for i in 1:3
    p = vline!(times, [5.0], color = :black, linestyle = :dot, label = "Drying", subplot = i)
    p = vline!(times, [10.0], color = :black, linestyle = :dash, label = "Intervention", subplot = i)
end
p
Plots.pdf("Drainage.pdf")

file = "data/CF_outline.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
active = Array(.!isnan.(Array(cf)))

sumabuns = Float64.(sum(abuns[numMoss+1:end, :, :, end], dims = 1))[1, :, :]
sumabuns[.!active] .= NaN
heatmap(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', size = (1000, 800), aspect_ratio = 1,
margin = 1.0*Plots.inch)
Plots.pdf("EndDrainage.pdf")
