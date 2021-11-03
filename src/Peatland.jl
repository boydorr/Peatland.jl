module Peatland

using Base: Int64
include("AbioticEnv.jl")
export peatlandAE

include("Scenarios.jl")
export Invasive, WindDispersal

include("AgeStructure.jl")
export Ageing

include("Dispersal.jl")
export PeatLookup

include("PeatSystem.jl")
export Ecosystem, update_peat_environment!

end
