module Peatland

using Base: Int64
include("AbioticEnv.jl")
export peatlandAE

include("Scenarios.jl")
export Invasive, WindDispersal

include("AgeStructure.jl")
export Ageing

end
