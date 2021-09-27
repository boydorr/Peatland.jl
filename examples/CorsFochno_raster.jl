using RCall
R"
library(rgdal)
library(rgeos)
library(raster)

cf <- readOGR('data/CorsFochno.shp')
r_cf <- raster(extent(cf), resolution = 50, crs = projection(cf))
r_cf <- rasterize(cf, r_cf)
plot(r_cf)
writeRaster(r_cf, 'data/CorsFochno.tif', overwrite=T)
"

