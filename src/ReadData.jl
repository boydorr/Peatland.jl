import EcoSISTEM.ClimatePref: UNITDICT
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using NetCDF

push!(UNITDICT, "hours since 1970-01-01T00:00:00Z" => Unitful.hr)
push!(UNITDICT, "seconds since 1961-01-01 00:00:00" => Unitful.s)
push!(UNITDICT, "kg m-2" => kg/m^2)

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



function readSoilCHESS(file::String, param::String; level = 1)
    x = reverse(ncread(file, "eastings"))
    y = ncread(file, "northings")
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1961years)
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
    array = ncread(file, param) * 1.0
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    
    return AxisArray(array * units, Axis{:x}(x * m), Axis{:y}(y * m),
                        Axis{:time}(time))
end
function readSoilCHESS(file::String, param::String, times::Vector{T}, xs::Vector{L}, ys::Vector{L}; level = 1) where {T <: Unitful.Time, L <: Unitful.Length}
    x = reverse(ncread(file, "eastings")) * m
    y = ncread(file, "northings") * m
    select_xs = (x .>= xs[1]) .& (x .<= xs[end])
    select_ys = (y .>= ys[1]) .& (y .<= ys[end])
    selects = findall(select_xs .& select_ys)
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1961years)
    select_times = findall((time .>= times[1]) .& (time .<= times[end]))
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
   
    array = NetCDF.readvar(NetCDF.open(file), param, start=[minimum(selects), level, minimum(select_times)],count = [length(selects), level, length(select_times)])
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    array = array * 1.0

    # reshape array and select level
    minX = minimum(x[selects]); maxX = maximum(x[selects]); minY = minimum(y[selects]); maxY = maximum(y[selects])
    array_new = zeros(Float64, length(minX:1000m:maxX), length(minY:1000m:maxY), length(time))
    array_new = AxisArray(array_new * units, Axis{:x}(minX:1000m:maxX), Axis{:y}(minY:1000m:maxY), Axis{:time}(time))
    newxs = x[selects]; newys = y[selects]
    for i in eachindex(newxs)
        for j in eachindex(time)
            east = newxs[i]; north = newys[i] 
            array_new[east .. east, north .. north, time[j] .. time[j]] .= array[i, 1, j] * units
        end
    end
    
    return array_new
end