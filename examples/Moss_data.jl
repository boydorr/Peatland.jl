using CSV
using DataFrames
using JLD2

moss_spp = CSV.File("data/Bryoatt_updated_2017/Bryoatt_updated_2017.csv", normalizenames=true)
moss_spp = DataFrame(moss_spp)
moss_spp = filter(x -> !ismissing(x.Taxon_name), moss_spp)

peat_spp = CSV.File("data/Peatland_species.csv", normalizenames=true)
peat_spp = DataFrame(peat_spp)

# Filter for species in the compiled peatland list
peat_spp_names = rstrip.(peat_spp[!, :Species_name_or_special_variable])
moss_spp = filter(row -> row.Taxon_name ∈ peat_spp_names, moss_spp)
# Filter for species found in Wales
moss_spp = filter(row -> !ismissing(row.W), moss_spp)
# Filter for moss
moss_spp = filter(row -> row.ML == "M", moss_spp)
# Filter for species that grow on peat
moss_spp = filter(row -> !ismissing(row.PT), moss_spp)

# Check if mosses are found in area according to NBN
cf_bryo = CSV.File("data/records-2022-09-01.csv", normalizenames=true)
cf_bryo = DataFrame(cf_bryo)
cf_spp = cf_bryo[!, :Species_Name]
moss_spp = filter(row -> row.Taxon_name ∈ cf_spp, moss_spp)
JLD2.save("data/Peat_30_moss.jld2", "moss_spp", moss_spp)

plant_spp = CSV.File("data/PLANTATT_19_Nov_08.csv", normalizenames=true)
plant_spp = DataFrame(plant_spp)
plant_spp = filter(x -> !ismissing(x.Taxon_name), plant_spp)
plant_spp = filter(x -> x.Taxon_name ∈ peat_spp.Species, plant_spp)


