module Peatland

using Base: Int64
include("AbioticEnv.jl")
export peatlandAE

include("Scenarios.jl")
export Invasive, WindDispersal, WaterUse, WaterFlux, Dry, 
Rewet, LateralFlow, changerate!

include("AgeStructure.jl")
export Ageing

include("Dispersal.jl")
export PeatLookup

include("PeatSystem.jl")
export PeatSystem, update_peat_environment!

include("Traits.jl")
export soiltrait, soilmatch

end
