using SimpleSDMLayers
using Plots
using Statistics

minLat = 289540      #lower latitude -> south
maxLat = 293461      #higher latitude -> north
minLong = 261527     #lower longitude -> west
maxLong = 265550     #higher longitude -> east

dir = "data\\2m_res_SN69_dsm"
fileList = readdir(dir)

function overlap(map)
    lat = extrema(latitudes(map))
    long = extrema(longitudes(map))

    inlat = lat[2] < maxLat  && lat[1] > minLat
    inlong = long[2] < maxLong && long[1] > minLong

    return inlat && inlong
end

mapList = SimpleSDMPredictor{Float64}[]
for file in eachindex(fileList)
    map = SimpleSDMLayers.ascii(joinpath(dir, fileList[file]))
    if overlap(map)
        map = coarsen(map, mean, (5, 5))
        push!(mapList, map)
    end
end

corsFochno = mosaic(mean, mapList)
corsFochno.grid[isnothing.(corsFochno.grid)] .= 0
corsFochno.grid[corsFochno.grid .< 0] .= 0
corsFochno.grid .= log.(corsFochno.grid)
plot(corsFochno)

using Shapefile
shp = Shapefile.shapes(Shapefile.Table("data/CorsFochno.shp"))
plot!(shp, fillcolor = false, linecolor = :white)

