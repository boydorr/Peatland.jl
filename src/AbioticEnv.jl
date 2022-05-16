using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
import EcoSISTEM: GridAbioticEnv, ContinuousHab, AbstractEcosystem, 
HabitatUpdate, cancel, checkbud, matchdict, NoChange


function Drying(eco::E, hab::ContinuousHab, timestep::Unitful.Time,
  thresh::Unitful.Volume = 0.0m^3) where E <: AbstractEcosystem
  val = hab.change.rate
  v = uconvert(m^3/unit(timestep), val)
  hab.matrix .+= v * timestep
  hab.matrix[hab.matrix .< thresh] .= thresh
end

function peatland_habitat(val::Unitful.Quantity, size::Unitful.Length, dim::Tuple{Int64, Int64})
    H = fill(val, dim)
    func = Drying
    rate = 0.0 * unit(val)/s
    habitatupdate = HabitatUpdate(func, rate, typeof(dimension(val)))
    return ContinuousHab(H, size, habitatupdate)
end

function elevation_habitat(ele::AxisArray, size::Unitful.Length)
  habitatupdate = HabitatUpdate(NoChange, 0.0/month)
  return ContinuousHab(ele.data, size, habitatupdate)
end

"""
     peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
         dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
         active::Array{Bool, 2})
"""
function peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
  active::Array{Bool, 2})
  if typeof(val) <: Unitful.Area
      val = uconvert(m^3, val)
  end
  area = uconvert(km^2, area)
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab = peatland_habitat(val, gridsquaresize, dimension)
  B = cancel(maxbud, area)
  bud = zeros(typeof(B), dimension)
  fill!(bud, B/(dimension[1]*dimension[2]))
  checkbud(B) || error("Unrecognised unit in budget")
  budtype = matchdict[unit(B)]
  return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end


function peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

  active = fill(true, dimension)
  peatlandAE(val, dimension, maxbud, area, active)
end

function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  maxbud::Unitful.Quantity{Float64})
  active = Array(.!isnan.(Array(ele)))
  if typeof(val) <: Unitful.Area
    val = uconvert(m^3, val)
  end
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab1 = peatland_habitat(val, gridsquaresize, dimension)
  hab2 = elevation_habitat(ele, gridsquaresize)
  hab = HabitatCollection2(hab1, hab2)

  B = cancel(maxbud, area)
  bud = zeros(typeof(B), dimension)
  fill!(bud, B/(dimension[1]*dimension[2]))
  checkbud(B) || error("Unrecognised unit in budget")
  budtype = matchdict[unit(B)]
  return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end



function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  maxbud::Unitful.Quantity{Float64}, active::Array{Bool, 2})
  if typeof(val) <: Unitful.Area
    val = uconvert(m^3, val)
  end
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab1 = peatland_habitat(val, gridsquaresize, dimension)
  hab2 = elevation_habitat(ele, gridsquaresize)
  hab = HabitatCollection2(hab1, hab2)

  B = cancel(maxbud, area)
  bud = zeros(typeof(B), dimension)
  fill!(bud, B/(dimension[1]*dimension[2]))
  checkbud(B) || error("Unrecognised unit in budget")
  budtype = matchdict[unit(B)]
  return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end


