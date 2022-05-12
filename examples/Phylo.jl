using PhyloNetworks
using PhyloPlots
using CSV
using DataFrames
using GLM
using JLD2

peat_spp = CSV.File("data/Peatland30spp.csv", normalizenames=true)
peat_spp = DataFrame(peat_spp)
peat_spp[!, :Plant_height_m_]

tree = readTopology("data/Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")

spp_names = peat_spp[!, :Species]
cross_species = spp_names ∩ tip_names
missing_species = setdiff(tip_names, cross_species)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end
peat_spp[!, :tipNames] = join.(split.(peat_spp[!, :Species], " "), "_")
peat_spp = filter(:Species => s -> s ∈ cross_species, peat_spp)

f = @formula(Plant_height_m_ ~ 1)
fitPagel = phylolm(f, peat_spp, tree, model="lambda")
lambda_estim(fitPagel)
recon = ancestralStateReconstruction(peat_spp[!, [:tipNames, :Plant_height_m_]], tree)

expect = expectations(recon)
tipName = [n.name for n in tree.leaf]
tipNumber = [n.number for n in tree.leaf]
tipOrder = indexin(peat_spp[!, :tipNames], tipName)
peat_spp[!, :tipNums] = tipNumber[tipOrder]
peat_spp[!, :Plant_height] = expect[indexin(peat_spp[!, :tipNums], expect[!, :nodeNumber]), :condExpectation]

peat_spp[!, :Leaf_width] = mean([peat_spp[!, :Leaf_min_width_cm] peat_spp[!, :Leaf_max_width_cm]], dims = 2)[:, 1]
peat_spp[!, :Leaf_length] = mean([peat_spp[!, :Leaf_min_length_cm] peat_spp[!, :Leaf_max_length_cm]], dims = 2)[:, 1]
f = @formula(Leaf_width ~ 1)
fitPagel = phylolm(f, peat_spp, tree, model="lambda")
lambda_estim(fitPagel)
recon = ancestralStateReconstruction(peat_spp[!, [:tipNames, :Leaf_width]], tree)
expect = expectations(recon)
peat_spp[!, :Leaf_width] = expect[indexin(peat_spp[!, :tipNums], expect[!, :nodeNumber]), :condExpectation]

f = @formula(Leaf_length ~ 1)
fitPagel = phylolm(f, peat_spp, tree, model="lambda")
lambda_estim(fitPagel)
recon = ancestralStateReconstruction(peat_spp[!, [:tipNames, :Leaf_length]], tree)
expect = expectations(recon)
peat_spp[!, :Leaf_length] = expect[indexin(peat_spp[!, :tipNums], expect[!, :nodeNumber]), :condExpectation]

f = @formula(Root_depth_max ~ 1)
fitdiscrete(tree, :ERSM, peat_spp[!, :Species], peat_spp[!, [:Species, :Root_depth_max]])
fitPagel = phylolm(f, peat_spp, tree, model="lambda")
lambda_estim(fitPagel)
recon = ancestralStateReconstruction(peat_spp[!, [:tipNames, :Root_depth_max]], tree)
expect = expectations(recon)
peat_spp[!, :Root_depth_max] = expect[indexin(peat_spp[!, :tipNums], expect[!, :nodeNumber]), :condExpectation]

PhyloPlots.plot(tree, :R, nodeLabel = expectationsPlot(recon))
JLD2.save("data/Peat_30_spp.jld2", "peat_spp", peat_spp)

