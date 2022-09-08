using EcoSISTEM
import EcoSISTEM: run_rule!

function run_rule!(eco::Ecosystem, rule::S, timestep::Unitful.Time) where S <: AbstractPeatState
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

function run_rule!(eco::Ecosystem, rule::P, timestep::Unitful.Time) where P <: AbstractPeatPlace
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

function run_rule!(eco::Ecosystem, rule::S, timestep::Unitful.Time) where S <: AbstractPeatSetUp
    if typeof(rule) == UpdateEnergy
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == WaterFlux
        _run_rule!(eco, rule, timestep)
    else
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::W, timestep::Unitful.Time) where W <: AbstractPeatWindDown
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
