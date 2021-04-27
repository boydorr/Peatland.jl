using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
import EcoSISTEM: GridAbioticEnv, ContinuousHab, AbstractEcosystem, HabitatUpdate, cancel, checkbud, matchdict

function Drying(eco::E, hab::ContinuousHab, timestep::Unitful.Time) where E <: AbstractEcosystem
  val = hab.change.rate
  v = uconvert(m^3/unit(timestep), val)
  hab.matrix .+= v * timestep
  hab.matrix[hab.matrix .< 0.0m^3] .= 0.0m^3
end

function peatland_habitat(val::Unitful.Quantity, size::Unitful.Length, dim::Tuple{Int64, Int64})
    H = fill(val, dim)
    func = Drying
    rate = 0.0 * unit(val)/s
    habitatupdate = HabitatUpdate(func, rate, typeof(dimension(val)))
    return ContinuousHab(H, size, habitatupdate)
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
