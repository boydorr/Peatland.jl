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
function runDrying(burnin::Unitful.Time, times1::Unitful.Time, times2::Unitful.Time, timestep::Unitful.Time, record_interval::Unitful.Time; save = false, save_path = pwd())
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
    bog_transitions = TransitionList(true)
    ditch_transitions = TransitionList(true)
    addtransition!(bog_transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(ditch_transitions, UpdateEnergy(update_energy_usage!))
    addtransition!(bog_transitions, UpdateEnvironment(update_peat_environment!))
    addtransition!(ditch_transitions, UpdateEnvironment(update_peat_environment!))
    active_squares = findall(active[1:end])
    peat_squares = findall(soil .== 11)
    peat_locs = [convert_coords(d[1], d[2], size(cf, 1)) for d in peat_squares]
    # Add in species based transitions (only for active squares)
    for spp in eachindex(sppl.species.names) 
        for loc in active_squares
            addtransition!(bog_transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(ditch_transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
            addtransition!(bog_transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(ditch_transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
            addtransition!(bog_transitions, SeedDisperse(spp, loc))
            addtransition!(ditch_transitions, SeedDisperse(spp, loc))
            # if spp > numMoss
                addtransition!(bog_transitions, Invasive(spp, loc, 1.0/30days))
                addtransition!(ditch_transitions, Invasive(spp, loc, 1.0/30days))
            # end
        end
    end
    # Water needs to be used everywhere (with a background rate for where we aren't modelling plants)
    for spp in eachindex(sppl.species.names) 
        for loc in eachindex(abenv.habitat.h1.matrix)
            addtransition!(bog_transitions, WaterUse(spp, loc, 1.0, 0.02, 0.9/100.0mm))
            addtransition!(ditch_transitions, WaterUse(spp, loc, 1.0, 0.02, 0.9/100.0mm))
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
        addtransition!(bog_transitions, LateralFlow(loc, κ, ν))
        addtransition!(bog_transitions, WaterFlux(loc, fmax, kₛ, W0, k2))
    end

    # Create ecosystem
    eco = PeatSystem(sppl, abenv, rel, transitions = bog_transitions)

    # Run simulation
    # Simulation Parameters
    lensim = length(0years:record_interval:burnin)
    lensim1 = length(0years:record_interval:times1)
    lensim2 = length(0years:record_interval:times2)
    abuns1 = generate_storage(eco, lensim1, 1)[:, :, :, 1]
    abuns2 = generate_storage(eco, lensim2, 1)[:, :, :, 1]
    # Burnin
    abuns = generate_storage(eco, lensim, 1)[:, :, :, 1]
    @time simulate_record!(abuns, eco, burnin, record_interval, timestep, save = save, save_path = save_path, specialise = true)
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    # Transitions for drying bog only
    file = "data/Ditches_full.tif"
    ditches = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
    ditch_locations = findall(.!isnan.(ditches))
    ditch_locs = [convert_coords(d[1], d[2], size(active, 1)) for d in ditch_locations]
    file = "data/MainRivers.tif"
    rivers = readfile(file, 289000.0m, 293000.0m, 261000.0m, 266000.0m)
    river_locations = findall(.!isnan.(rivers))
    river_locs = [convert_coords(r[1], r[2], size(active, 1)) for r in river_locations]
    drains = [ditch_locs; river_locs ...]
    abenv.habitat.h1.matrix[drains] .= 0.0
    abenv.habitat.h2.matrix[drains] .= 0
    abenv.active[ditch_locations] .= false
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
        addtransition!(ditch_transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(ditch_transitions, Drainage(loc, drainage, κ, ν))
    end
    # Neighbour squares have their water pulled into the ditch
    for loc in neighbour_drains
        drainage = 0.0/month
        κ = 1 * 3.0m^2/month
        ν = 1 * 3.0m^2/month
        addtransition!(ditch_transitions, LateralFlow(loc, κ, ν))
        κ = 0.001m^2/month
        ν = 0.001m^2/month
        addtransition!(ditch_transitions, Drainage(loc, drainage, κ, ν))
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
        addtransition!(ditch_transitions, LateralFlow(loc, κ, ν))
        addtransition!(ditch_transitions, WaterFlux(loc, fmax, kₛ, W0, k2))
        # But also experience a certain amount of drying
        addtransition!(ditch_transitions, Dry(loc, 0.3, 1year))
    end

    eco.transitions = ditch_transitions
    @time simulate_record!(abuns1, eco, times1, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    abenv.habitat.h2.matrix[drains] .= ele[drains]
    abenv.active[drains] .= true
    eco.transitions = bog_transitions

    @time simulate_record!(abuns2, eco, times2, record_interval, timestep, save = save, save_path = save_path, specialise = true);
    println(mean(eco.abenv.habitat.h1.matrix[active]))

    abuns = cat(abuns, abuns1, abuns2,  dims = 3)
    abuns = reshape(abuns, (numSpecies, grid[1], grid[2], lensim + lensim1 + lensim2))
    
    return abuns
end

burnin = 70year; times1 = 30year; times2 = 30year; timestep = 1month; record_interval = 3months
abuns = runDrying(burnin, times1, times2, timestep, record_interval);
@save "/home/claireh/sdc/Peatland/Peatland_baseline_ditches_new.jld2" abuns=abuns[:, :, :, [50*4, 70*4, 90*4, end]]


times = ustrip.(uconvert.(year, 20years:record_interval:((burnin+times1+times2)+record_interval*3)))
plot(grid = false, label = "", layout = (3, 1), left_margin = 1.0*Plots.inch, size = (1000, 1200))
for i in mosses
    plot!(times,sum(abuns[i, :, :, 80:end], dims = (1,2))[1, 1, :,], grid = false, label = "", subplot = 1,
    title = "Mosses", colour = :grey, alpha = 0.1)
end
plot!(times, mean(sum(abuns[mosses, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 1,
    colour = :black)
vline!([burnin/year, (burnin+times1)/year], subplot = 1, color = :black, lty = :dot, label = "")
for i in shrubs
    plot!(times,sum(abuns[i, :, :, 80:end], dims = (1, 2))[1, 1, :,], grid = false, label = "", subplot = 2,
    title = "Shrubs, Herbs and Grasses", colour = :grey, alpha = 0.1)
end
plot!(times, mean(sum(abuns[shrubs, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 2,
     colour = :black)
vline!([burnin/year, (burnin+times1)/year], subplot = 2, color = :black, lty = :dot, label = "")
for i in trees
    plot!(times,sum(abuns[i, :, :, 80:end], dims = (1, 2))[1, 1, :,], grid = false, label = "", subplot = 3, 
    title = "Trees", colour = :grey, alpha = 0.1)
end
plot!(times, mean(sum(abuns[trees, :, :, 80:end], dims = (2, 3))[:, 1, 1, :], dims = 1)[1, :], grid = false, label = "", subplot = 3,
   colour = :black)
vline!([burnin/year, (burnin+times1)/year], subplot = 3, color = :black, lty = :dot, label = "")
Plots.pdf("Abuns_baseline_ditches_new.pdf")


@load "data/Peatland_baseline_ditches_new.jld2"
JLD2.@load("data/Peat_30_spp.jld2", peat_spp)
JLD2.@load("data/Peat_30_moss.jld2", moss_spp)
numMoss = nrow(moss_spp); numShrub = nrow(peat_spp); numSpecies = numMoss + numShrub
mosses = 1:numMoss
shrubs = findall((peat_spp[!, :Type] .== "Shrub") .| (peat_spp[!, :Type] .== "Grass") .| (peat_spp[!, :Type] .== "Herb")) .+ numMoss
trees = findall((peat_spp[!, :Type] .== "Tree")) .+ numMoss
file = "data/CF_outline.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
active = Array(.!isnan.(Array(cf)))

heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 600), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
sumabuns = Float64.(reshape(sum(abuns[mosses, :, :, 1], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A")
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, :, 2], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B")
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, :, 3], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C")
sumabuns = Float64.(reshape(sum(abuns[1:numMoss, :, :, end], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Mosses_baseline_new.pdf")

heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 100), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, 1], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A")
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, 2], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B")
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, 3], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C")
sumabuns = Float64.(reshape(sum(abuns[shrubs, :, :, end], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Shrubs_baseline_new.pdf")

heatmap(layout = (2,2), size = (1200, 1000), grid = false, aspect_ratio = 1, clim = (0, 2), titlelocation = :left,
left_margin = 10*Plots.mm, guidefont = "Helvetica Bold", guidefontsize = 16, titlefontsize = 18)
sumabuns = Float64.(reshape(sum(abuns[trees, :, :, 1], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 1, c = :viridis,
title = "A")
sumabuns = Float64.(reshape(sum(abuns[trees, :, :, 2], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 2, c = :viridis,
title = "B")
sumabuns = Float64.(reshape(sum(abuns[trees, :, :, 3], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 3, c = :viridis,
title = "C")
sumabuns = Float64.(reshape(sum(abuns[trees, :, :, end], dims = 1), size(active)))
sumabuns[.!active] .= NaN
heatmap!(261000.0:10:266000.0, 289000.0:10:293000.0, sumabuns', subplot = 4, c = :viridis,
title = "D")
Plots.pdf("plots/Trees_baseline_new.pdf")