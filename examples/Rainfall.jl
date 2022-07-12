using CSV
using DataFrames
using Dates
using Plots
using JLD2

rainfall = CSV.read("data/RainfallCF.csv", DataFrame)
missing_rows = ismissing.(rainfall[!, "AWS CF"])
rainfall[missing_rows, "AWS CF"] .= rainfall[missing_rows, "CF correction (88.3%)"]

rainfall[!, :DATE] .= Date.(rainfall[!, :DATE], "dd/mm/yyyy")
rainfall[!, :month] = Dates.month.(rainfall[!, :DATE])
rainfall[!, :year] = Dates.year.(rainfall[!, :DATE])
gdf = groupby(rainfall, [:year, :month])
monthly_rain = combine(gdf, "AWS CF" => (x -> sum(skipmissing(x))) => :sum)
plot(monthly_rain[!, :sum])

file = "data/CorsFochno.tif"
cf = readfile(file, 261659.3m, 265409.3m, 289536.7m, 292756.7m)
grid = size(cf)
prec_array = fill(0.0mm, grid[1], grid[2], nrow(monthly_rain))
for i in 1:nrow(monthly_rain)
    prec_array[:, :, i] .= fill(monthly_rain[i, :sum], grid) * mm
end

bud = WaterTimeBudget(prec_array, 1)
@save "data/RainfallBudget.jld2" bud