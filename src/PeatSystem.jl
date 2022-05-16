import EcoSISTEM: Ecosystem, AbstractTraitRelationship, 
emptygridlandscape, create_cache, AbstractCache, numrequirements,
_invalidatecaches!, invalidatecaches!, budgetupdate!, habitatupdate!, genlookups, getkernels, SpeciesLookup,
_getdimension
using Distributions
"""
  PeatCache

Cache houses an integer array of seed production made by all species in the Ecosystem,
`seedbank`, moves made by all dispersing seeds in a timestep, `netmigration`,
a matrix of current resource used in the Ecosystem, `totalE`,
and a Bool to say if these caches are `valid`.
"""
mutable struct PeatCache <: AbstractCache
  netmigration::Array{Float64, 2}
  seedbank::Array{Int64, 2}
  watermigration::Array{Unitful.Volume, 2}
  totalE::Matrix{Float64}
  valid::Bool
end

function create_peat_cache(abenv::A, sppl::SpeciesList{SpeciesTypes{TR, R, MO, T}},
  ml::GridLandscape) where {A <: AbstractAbiotic, TR, R, MO <: BirthOnlyMovement, T}
nm = zeros(Float64, size(ml.matrix))
sb = zeros(Int64, size(ml.matrix))
wm = zeros(typeof(1.0m^3), _getdimension(abenv.habitat))
totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
return PeatCache(nm, sb, wm, totalE, false)
end

function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
    rel::AbstractTraitRelationship, lookup::PeatLookup; transitions::Union{Nothing, TransitionList} = nothing) where {F<:Function, T, Req}
 
     # Create matrix landscape of zero abundances
   ml = emptygridlandscape(abenv, spplist)
   # Populate this matrix with species abundances
   popfun(ml, spplist, abenv, rel)
   cache = create_peat_cache(abenv, spplist, ml)

   return Ecosystem{typeof(ml), typeof(abenv), typeof(spplist), typeof(rel), typeof(lookup), typeof(cache)}(ml, spplist, abenv,
   missing, rel, lookup, cache, transitions)
 end

 function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
    rel::AbstractTraitRelationship, lookup::PeatLookup; transitions::Union{Nothing, TransitionList} = nothing)
    return Ecosystem(populate!, spplist, abenv, rel, lookup, transitions = transitions)
 end

 function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
  rel::AbstractTraitRelationship; transitions::Union{Nothing, TransitionList} = nothing, peatcache::Bool = true) where {F<:Function, T, Req}

   # Create matrix landscape of zero abundances
 ml = emptygridlandscape(abenv, spplist)
 # Populate this matrix with species abundances
 popfun(ml, spplist, abenv, rel)
 if peatcache
    cache = create_peat_cache(abenv, spplist, ml)
 end
 lookup = SpeciesLookup(collect(map(k -> EcoSISTEM.genlookups(abenv.habitat, k), getkernels(spplist.species.movement))))
 
 return Ecosystem{typeof(ml), typeof(abenv), typeof(spplist), typeof(rel), typeof(lookup), typeof(cache)}(ml, spplist, abenv,
 missing, rel, lookup, cache, transitions)
end

function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
  rel::AbstractTraitRelationship; transitions::Union{Nothing, TransitionList} = nothing, peatcache::Bool = true)
  return Ecosystem(populate!, spplist, abenv, rel, transitions = transitions, peatcache = peatcache)
end

 function _invalidatecaches!(eco::Ecosystem, cache::PeatCache)
  eco.ordinariness = missing
  eco.cache.netmigration .= 0
  eco.cache.seedbank .= 0
  eco.cache.watermigration .= 0m^3
  eco.cache.valid = false
end

 function update_peat_environment!(eco::Ecosystem, timestep::Unitful.Time)
  rng = eco.abundances.rngs[Threads.threadid()]
  eco.abundances.matrix .+= rand.(rng, Poisson.(eco.cache.netmigration))
  
  hab = getsoilwater(eco.abenv)
  hab .+= eco.cache.watermigration

  # Invalidate all caches for next update
  invalidatecaches!(eco)

  # Update environment - habitat and energy budgets
  habitatupdate!(eco, timestep)
  budgetupdate!(eco, timestep)
end
