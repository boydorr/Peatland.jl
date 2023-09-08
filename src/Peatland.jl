module Peatland

using Base: Int64
include("AbioticEnv.jl")
export peatlandAE

include("Scenarios.jl")
export Invasive, WindDispersal, WaterUse, WaterFlux, Dry, 
Rewet, LateralFlow, changerate!, Drainage

include("AgeStructure.jl")
export Ageing

include("Dispersal.jl")
export PeatLookup

include("PeatSystem.jl")
export PeatSystem, update_peat_environment!

include("Traits.jl")
export soiltrait, soilmatch

# include("Run.jl")
# export run_rule!

include("ReadData.jl")
export readCHESS, readSoilCHESS
end
