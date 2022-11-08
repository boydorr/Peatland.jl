using Revise
using CSV
using DataFrames
using Dates
using Plots
using JLD2
using Peatland
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM
using EcoSISTEM.Units
import EcoSISTEM.Units: month, year
using EcoSISTEM.ClimatePref
using Statistics

rainfall = CSV.read("data/RainfallCF.csv", DataFrame)
missing_rows = ismissing.(rainfall[!, "AWS CF"])
rainfall[missing_rows, "AWS CF"] .= rainfall[missing_rows, "CF correction (88.3%)"]

rainfall[!, :DATE] .= Date.(rainfall[!, :DATE], "dd/mm/yyyy")
rainfall[!, :month] = Dates.month.(rainfall[!, :DATE])
rainfall[!, :year] = Dates.year.(rainfall[!, :DATE])
gdf = groupby(rainfall, [:year, :month])
monthly_rain = combine(gdf, "AWS CF" => (x -> sum(skipmissing(x))) => :sum)
plot(monthly_rain[!, :sum])

file = "data/CF_outline.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
grid = size(cf)
prec_array = fill(0.0mm, grid[1], grid[2], nrow(monthly_rain))
for i in 1:nrow(monthly_rain)
    prec_array[:, :, i] .= fill(monthly_rain[i, :sum], grid) * mm
end

bud = WaterTimeBudget(prec_array, 1)
@save "data/RainfallBudget.jld2" bud

using NetCDF
file = "data/chess-scape_rcp85_bias-corrected_01_pr_uk_1km_monthly_19801201-20801130.nc"
param = "pr"
times = 2010years+1month:1Units.month:2080year
xs = 261000.0m:1m:266000.0m; ys = 289000.0m:1m:293000.0m
pr = readCHESS(file, param, collect(times), collect(xs), collect(ys))

# 1kg water per 1m^2 is 1mm thick
pr = uconvert.(mm, (ustrip.(pr) .* mm/s) .* 1month)
pr = mean(pr, dims = (1, 2))[1, 1, :]
plot!(1:843, ustrip.(pr))

prec_array2 = fill(0.0mm, grid[1], grid[2], length(pr))
for i in eachindex(pr)
    prec_array2[:, :, i] .= pr[i]
end
bud = WaterTimeBudget(prec_array, 1)
@save "data/RainfallBudget_future.jld2" bud

bud = WaterTimeBudget(prec_array2[:, :, 1:120], 1)
@save "data/RainfallBudget_burnin.jld2" bud