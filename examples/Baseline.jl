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
function runDrying(timestep::Unitful.Time; save = false, save_path = pwd())
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
    file = "data/LCM.tif"
    soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
    soil = Int.(soil)
    abenv = peatlandAE(tpi, ele, soil, bud, active)

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

    # Add in location based transitions
    for loc in active_squares
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

    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = transitions)

    # Run simulation
    # Simulation Parameters
    burnin = 30year; times1 = 10year; times2 = 10year; 
    record_interval = 3months
    lensim = length(0years:record_interval:burnin)
    lensim1 = length(0years:record_interval:times1)
    lensim2 = length(0years:record_interval:times2)
    abuns1 = generate_storage(eco, lensim1, 1)[:, :, :, 1]
    abuns2 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate_record!(abuns, eco, burnin, record_interval, timestep, save = save, save_path = save_path, specialise = true)
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    for loc in active_squares
        addtransition!(transitions, Peatland.Dry(loc, 0.3, 1year))
    end
    @time simulate_record!(abuns1, eco, times1, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    for loc in active_squares
        addtransition!(transitions, Peatland.Rewet(loc, 0.3, 1year, 1.0))
    end
    @time simulate_record!(abuns2, eco, times2, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    abuns = cat(abuns, abuns1, abuns2,  dims = 3)
    abuns = reshape(abuns, (numSpecies, grid[1], grid[2], lensim + lensim1 + lensim2))
    
    return abuns
end

abuns = runDrying(1month);
@save "/home/claireh/sdc/Peatland/Peatland_baseline_new.jld2" abuns=abuns[:, :, :, [12,end]]


plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(0:3/12:30.75,sum(abuns[i, :, :, 80:end], dims = (1,2))[1, 1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
annotate!([-10], [7e6], ["A"])
plot!(0:3/12:30.75, mean(sum(abuns[mosses, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
vline!([10, 20], subplot = 1, color = :black, lty = :dot, label = "")
for i in shrubs
    plot!(0:3/12:30.75,sum(abuns[i, :, :, 80:end], dims = (1, 2))[1, 1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:30.75, mean(sum(abuns[shrubs, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
vline!([10, 20], subplot = 2, color = :black, lty = :dot, label = "")
for i in trees
    plot!(0:3/12:30.75,sum(abuns[i, :, :, 80:end], dims = (1, 2))[1, 1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(0:3/12:30.75, mean(sum(abuns[trees, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
vline!([10, 20], subplot = 3, color = :black, lty = :dot, label = "")
Plots.pdf("Abuns_baseline_new.pdf")


@load "data/Peatland_baseline_new.jld2"
JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
JLD2.@load("data/Peat_30_moss.jld2", moss_spp)
numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub
mosses = 1:numMoss
shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss
file = "data/CF_outline.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
active = Array(.!isnan.(Array(cf)))

heatmap(layout = (1,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 10), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, 1], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A")
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, 2], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B")
Plots.pdf("plots/Shrubs.pdf")
