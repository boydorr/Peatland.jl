import EcoSISTEM.ClimatePref: UNITDICT
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using NetCDF

push!(UNITDICT, "hours since 1970-01-01T00:00:00Z" => Unitful.hr)
const NEWUNITDICT = UNITDICT
function readCHESS(file::String, param::String)
    x = reverse(ncread(file, "x"))
    y = ncread(file, "y")
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1970years)
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
    array = ncread(file, param) * 1.0
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    
    return AxisArray(array * units, Axis{:x}(x * m), Axis{:y}(y * m),
                        Axis{:time}(time))
end
function readCHESS(file::String, param::String, times::Vector{T}, xs::Vector{L}, ys::Vector{L}) where {T <: Unitful.Time, L <: Unitful.Length}
    x = reverse(ncread(file, "x")) * m
    y = ncread(file, "y") * m
    select_xs = findall((x .>= xs[1]) .& (x .<= xs[end]))
    select_ys = findall((y .>= ys[1]) .& (y .<= ys[end]))
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1970years)
    select_times = findall((time .>= times[1]) .& (time .<= times[end]))
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
   
    array = NetCDF.readvar(NetCDF.open(file), param, start=[minimum(select_xs), minimum(select_ys), minimum(select_times)],count = [length(select_xs),length(select_ys), length(select_times)])
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    array = array * 1.0
    
    return AxisArray(array * units, Axis{:x}(x[select_xs]), Axis{:y}(y[select_ys]),
                        Axis{:time}(time[select_times]))
end