import EcoSISTEM: ContinuousTrait

"""
    GaussTrait{C <: Number} <: ContinuousTrait{C}

Trait type that holds Gaussian mean and variance trait information for each species, of any number type `C`.
"""
function GaussTrait(mean::Array{C, 1}, var::Array{C, 1}) where C  <: Unitful.Area
    meanvol = uconvert.(m^3, mean)
    varvol = ustrip.(var) .* m^3
    return GaussTrait{typeof(1.0m^3)}(meanvol, varvol)
end

using EcoSISTEM
import EcoSISTEM: AbstractTraitRelationship, AbstractTraits, iscontinuous, _traitfun, getpref, combineTR
import Base.eltype

mutable struct soilmatch{TR} <: AbstractTraitRelationship{TR}
end
function (::soilmatch{TR})(niche::TR, pref::Vector{TR}) where TR
  if niche in pref
    return 1.0
  else
    return 0.0
  end
end
iscontinuous(tr::soilmatch{TR}) where TR = false
function eltype(tr::soilmatch{TR}) where TR
    return TR
end


mutable struct soiltrait{D <: Number} <: AbstractTraits{D}
  vals::Array{Array{D, 1}, 1}
end

iscontinuous(trait::soiltrait{D}) where D = false
function eltype(trait::soiltrait{D}) where D
    return D
end

function soiltrait(vals::Array{Array{D, 1}, 1}) where D  <: Integer
    return soiltrait{typeof(1)}(vals)
end

function _traitfun(hab::DiscreteHab, trts::soiltrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        vals = getpref(trts, spp)
    return rel(h, vals)
end

function getpref(traits::soiltrait, spp::Int64)
  return traits.vals[spp]
end