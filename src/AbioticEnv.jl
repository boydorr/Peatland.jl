using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
import EcoSISTEM: GridAbioticEnv, ContinuousHab, AbstractEcosystem, 
HabitatUpdate, cancel, checkbud, matchdict, NoChange, AbstractHabitat, AbstractBudget


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

function soil_habitat(soil::Array{Int, 2}, size::Unitful.Length)
  habitatupdate = HabitatUpdate(NoChange, 0.0/month)
  return DiscreteHab(soil, size, habitatupdate)
end

"""
     peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
         dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
         active::Array{Bool, 2})
"""
function peatlandAE(hab::H, active::Array{Bool, 2}, bud::B) where {H <: AbstractHabitat, B <: AbstractBudget}
  return GridAbioticEnv{H, B}(hab, active, bud)
end

function _habitat_builder(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, area::Unitful.Area{Float64})
  if typeof(val) <: Unitful.Area
      val = uconvert(m^3, val)
  end
  area = uconvert(km^2, area)
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  return peatland_habitat(val, gridsquaresize, dimension)
end

function _habitat_builder(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}})
  if typeof(val) <: Unitful.Area
    val = uconvert(m^3, val)
  end
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab1 = peatland_habitat(val, gridsquaresize, dimension)
  hab2 = elevation_habitat(ele, gridsquaresize)

  return HabitatCollection2(hab1, hab2)
end

function _habitat_builder(ele::AxisArray, soil::Array{Int64}, val::Union{Float64, Unitful.Quantity{Float64}})
  if typeof(val) <: Unitful.Area
    val = uconvert(m^3, val)
  end
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab1 = peatland_habitat(val, gridsquaresize, dimension)
  hab2 = elevation_habitat(ele, gridsquaresize)
  hab3 = soil_habitat(soil, gridsquaresize)

  return HabitatCollection3(hab1, hab2, hab3)
end

function _habitat_builder(tpi::AxisArray, ele::AxisArray, soil::Array{Int64})
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

  hab1 = elevation_habitat(tpi, gridsquaresize)
  hab2 = elevation_habitat(ele, gridsquaresize)
  hab3 = soil_habitat(soil, gridsquaresize)

  return HabitatCollection3(hab1, hab2, hab3)
end

function _budget_builder(dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})
  B = cancel(maxbud, area)
  bud = zeros(typeof(B), dimension)
  fill!(bud, B/(dimension[1]*dimension[2]))
  checkbud(B) || error("Unrecognised unit in budget")
  budtype = matchdict[unit(B)]
  return budtype(bud)
end

function _budget_builder(ele::AxisArray, maxbud::Unitful.Quantity{Float64})
  dimension = size(ele)
  area = uconvert(km^2, step(ele.axes[1]) * step(ele.axes[2]) * prod(dimension))
  B = cancel(maxbud, area)
  bud = zeros(typeof(B), dimension)
  fill!(bud, B/(dimension[1]*dimension[2]))
  checkbud(B) || error("Unrecognised unit in budget")
  budtype = matchdict[unit(B)]
  return budtype(bud)
end

function peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
  active::Array{Bool, 2})
  hab = _habitat_builder(val, dimension, area)
  bud = _budget_builder(dimension, maxbud, area)
  return peatlandAE(hab, active, bud)
end


function peatlandAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})
  active = fill(true, dimension)
  return peatlandAE(val, dimension, maxbud, area, active)
end

function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  maxbud::Unitful.Quantity{Float64})
  active = Array(.!isnan.(Array(ele)))
  hab = _habitat_builder(ele, val)
  bud = _budget_builder(ele, maxbud)
  return peatlandAE(hab, active, bud)
end

function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  maxbud::Unitful.Quantity{Float64}, active::Array{Bool, 2})
  hab = _habitat_builder(ele, val)
  bud = _budget_builder(ele, maxbud)
  return peatlandAE(hab, active, bud)
end


function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  bud::B) where B <: AbstractBudget
  active = Array(.!isnan.(Array(ele)))
  hab = _habitat_builder(ele, val)
  return peatlandAE(hab, active, bud)
end

function peatlandAE(ele::AxisArray, val::Union{Float64, Unitful.Quantity{Float64}}, 
  bud::B, active::Array{Bool, 2}) where B <: AbstractBudget
  hab = _habitat_builder(ele, val)
  return peatlandAE(hab, active, bud)
end

function peatlandAE(ele::AxisArray, soil::Array{Int64, 2}, val::Union{Float64, Unitful.Quantity{Float64}}, 
  bud::B) where B <: AbstractBudget
  active = Array(.!isnan.(Array(ele)))
  hab = _habitat_builder(ele, soil, val)
  return peatlandAE(hab, active, bud)
end

function peatlandAE(ele::AxisArray, soil::Array{Int64, 2}, val::Union{Float64, Unitful.Quantity{Float64}}, 
  bud::B, active::Array{Bool, 2}) where B <: AbstractBudget
  hab = _habitat_builder(ele, soil, val)
  return peatlandAE(hab, active, bud)
end


function peatlandAE(tpi::AxisArray, ele::AxisArray, soil::Array{Int64, 2}, bud::B) where B <: AbstractBudget
  active = Array(.!isnan.(Array(ele)))
  hab = _habitat_builder(tpi, ele, soil)
  return peatlandAE(hab, active, bud)
end

function peatlandAE(tpi::AxisArray, ele::AxisArray, soil::Array{Int64, 2}, bud::B, active::Array{Bool, 2}) where B <: AbstractBudget
  hab = _habitat_builder(tpi, ele, soil)
  return peatlandAE(hab, active, bud)
end
