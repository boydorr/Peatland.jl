import Simulation: ContinuousTrait

"""
    GaussTrait{C <: Number} <: ContinuousTrait{C}

Trait type that holds Gaussian mean and variance trait information for each species, of any number type `C`.
"""
function GaussTrait(mean::Array{C, 1}, var::Array{C, 1}) where C  <: Unitful.Area
    meanvol = uconvert.(m^3, mean)
    varvol = ustrip.(var) .* m^3
    return GaussTrait{typeof(1.0m^3)}(meanvol, varvol)
end
