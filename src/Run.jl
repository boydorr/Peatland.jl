using EcoSISTEM
import EcoSISTEM: run_rule!, AbstractStateTransition, AbstractPlaceTransition, AbstractSetUp, AbstractWindDown


# macro create_code(rule_types)

#     first_line = ["if typeof(rule) == $(rule_types[1]) _run_rule!(eco, rule, timestep)"]
#     for i in rule_types[2:end]
#         push!(first_line, "elseif typeof(rule) == $i _run_rule!(eco, rule, timestep)")
#     end
#     last_line = "else run_rule!(eco, rule, timestep) end"
#     body = Meta.parse(prod(first_line) * last_line)
#     return body
# end

function run_rule!(eco::Ecosystem, rule::AbstractStateTransition, timestep::Unitful.Time)
    if typeof(rule) == BirthProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == GenerateSeed
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == DeathProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Invasive
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == WaterUse
        _run_rule!(eco, rule, timestep)
    else 
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::AbstractPlaceTransition, timestep::Unitful.Time)
    if typeof(rule) == AllDisperse
        _run_rule!(eco, rule)
    elseif typeof(rule) == SeedDisperse
        _run_rule!(eco, rule)
    elseif typeof(rule) == WindDispersal
        _run_rule!(eco, rule)
    else
        _run_rule!(eco, rule)
    end
end

function run_rule!(eco::Ecosystem, rule::AbstractSetUp, timestep::Unitful.Time)
    if typeof(rule) == UpdateEnergy
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == WaterFlux
        _run_rule!(eco, rule, timestep)
    else
        _run_rule!(eco, rule, timestep)
    end
end
function run_rule!(eco::Ecosystem, rule::AbstractWindDown, timestep::Unitful.Time)
    if typeof(rule) == UpdateEnvironment
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Dry
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Rewet
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == LateralFlow
        _run_rule!(eco, rule, timestep)    
    else
        _run_rule!(eco, rule, timestep)
    end
end
