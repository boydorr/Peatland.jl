import EcoSISTEM: Ecosystem, AbstractTraitRelationship, emptygridlandscape, create_cache

function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
    rel::AbstractTraitRelationship, lookup::PeatLookup; transitions::Union{Nothing, TransitionList} = nothing) where {F<:Function, T, Req}
 
     # Create matrix landscape of zero abundances
   ml = emptygridlandscape(abenv, spplist)
   # Populate this matrix with species abundances
   popfun(ml, spplist, abenv, rel)
   cache = create_cache(spplist, ml)
   return Ecosystem{typeof(ml), typeof(abenv), typeof(spplist), typeof(rel), typeof(lookup), typeof(cache)}(ml, spplist, abenv,
   missing, rel, lookup, cache, transitions)
 end

 function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
    rel::AbstractTraitRelationship, lookup::PeatLookup; transitions::Union{Nothing, TransitionList} = nothing)
    return Ecosystem(populate!, spplist, abenv, rel, lookup, transitions = transitions)
 end