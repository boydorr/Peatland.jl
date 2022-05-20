using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols

import EcoSISTEM: getprob, getspecies, getlocation

@testset "Basic scenarios" begin
    spp = 1; loc = 1; prob = 1.0/day
    
    rule = Invasive(spp, loc, prob)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
    @test getprob(rule) == prob

    rule = WaterFlux(loc, prob)
    @test getlocation(rule) == loc
    @test rule.prob == prob

    len = 1year; time = 1month; flprob = 0.1
    rule = Dry(loc, flprob, len, time)
    @test getlocation(rule) == loc
    @test rule.prob == flprob
    rule = Dry(loc, flprob, len)
    @test rule.time == 1month

    maxvol = 10.0m^3
    rule = Rewet(loc, flprob, len, time, maxvol)
    @test getlocation(rule) == loc
    @test rule.prob == flprob
    rule = Rewet(loc, flprob, len, maxvol)
    @test rule.time == 1month

    frac = 0.1
    rule = WaterUse(spp, loc, frac)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
    @test rule.soil_moisture_frac == frac

    rule = WindDispersal(spp, loc, 1.0m, 0.3, 1.0m, 1.0m/s, 1.0m/s)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
end
