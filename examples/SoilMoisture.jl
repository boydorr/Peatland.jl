using NetCDF
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Peatland
using Statistics

times = collect(2015years:1month:2016years)
eastings = unique(reverse(ncread("data/CHESSLandMonHydEn2015.nc", "eastings"))) .* m
northings = unique(ncread("data/CHESSLandMonHydEn2015.nc", "northings")) .* m
xs = collect(261000.0m:1000m:266000.0m); ys = collect(289000.0m:1000m:293000.0m)
soilm = readSoilCHESS("data/CHESSLandMonHydEn2015.nc", "smcl", times, xs, ys)
mean(soilm) .* 10m^2 / 100kg

## Approx. 3.77 m^3

gmax = 0.05mm/(g*day*m^2)
gmax * 10m^2 * 1000g/m