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

gr()
file = "data/LCM.tif"
colours = ["red", "darkgreen", "brown", "green", "lightgreen", "green3", "yellow4", "yellow",
            "purple", "pink", "turquoise3", "lavender", "navy", "blue", "gold", "gold", "lightyellow","lightyellow", "lightblue", "black", "grey"]
labels = ["Broadleaved woodland", "Coniferous woodland", "Arable",
"Improved grassland", "Neutral grassland", "Calcareous grassland", "Acid grassland",
"Fen, Marsh, Swamp", "Heather",
"Heather grassland", "Bog",
"Inland Rock", "Saltwater",
"Freshwater", "Supra-littoral rock",
"Supra-littoral sediment", "Littoral rock",
"Littoral sediment", "Saltmarsh", "Urban",
"Suburban"]
soil = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
heatmap(soil', c = colours)
heatmap!(ditches', c = :white)

file = "data/CF_elevation.tif"
ele = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
heatmap(log.(1 .+ ele)',  cbar = true)
file = "data/Ditches2.tif"
ditches = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
heatmap!(ditches', c = :white)
