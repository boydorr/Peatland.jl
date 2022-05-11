using SimpleSDMLayers
using Plots
using Statistics
using RCall
using Shapefile
R"
library(rgdal)
library(rgeos)
library(raster)

cf <- readOGR('../../data/CorsFochno.shp')
ext <- matrix(extent(cf))
print(projection(cf))
"
@rget ext
minLat = ext[3]     #lower latitude -> south
maxLat = ext[4]    #higher latitude -> north
minLong = ext[1]   #lower longitude -> west
maxLong = ext[2] #higher longitude -> east

dir = "../../data/2m_res_dsm"
fileList = readdir(dir)

function overlap(map::SimpleSDMPredictor, latMin::Float64, latMax::Float64, 
    lonMin::Float64, lonMax::Float64)
    lat = extrema(latitudes(map))
    long = extrema(longitudes(map))

    inlat = (lat[1] <= latMax)  && (lat[2] >= latMin)
    inlong = (long[1] <= lonMax) && (long[2] >= lonMin)

    return inlat && inlong
end

mapList = SimpleSDMPredictor{Float64}[]
for file in eachindex(fileList)
    map = SimpleSDMLayers.ascii(joinpath(dir, fileList[file]))
    if overlap(map, minLat, maxLat, minLong, maxLong)
        map = coarsen(map, mean, (5, 5))
        push!(mapList, map)
    end
end

corsFochno = mosaic(mean, mapList)
corsFochno.grid[isnothing.(corsFochno.grid)] .= 0
plot(corsFochno, seriescolor = cgrad(:delta, [0, 1/200, 1]), clim = (0,200))

cf = Shapefile.shapes(Shapefile.Table("../../data/CorsFochno.shp"))
cf_poly = SimpleSDMLayers._format_polygon(GeoInterface.coordinates(cf[1])[1])
plot!(cf, fillcolor = false, linecolor = :black)
corsFochno_clip = mask(cf_poly, corsFochno)

plot(corsFochno_clip, seriescolor = :deep, grid = false)
plot!(cf, fillcolor = false, linecolor = :black)
Plots.pdf("../../plots/CorsFochno_elevation.pdf")

geotiff("../../data/CF_elevation.tif", corsFochno_clip)