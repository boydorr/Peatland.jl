using EcoSISTEM
using EcoSISTEM.Units
using Peatland
using Unitful, Unitful.DefaultSymbols

@testset "Basic scenarios" begin
    spp = 1; loc = 1; prob = 1.0/day
    
    rule = Invasive(spp, loc, prob)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
    @test getprob(rule) == prob

    rule = WaterFlux(loc, prob)
    @test getlocation(rule) == loc
    @test getprob(rule) == prob

    len = 1year; time = 1month
    rule = Dry(loc, prob, len, time)
    @test getlocation(rule) == loc
    @test getprob(rule) == prob

    len = 1year; time = 1month
    rule = Rewet(loc, prob, len, time)
    @test getlocation(rule) == loc
    @test getprob(rule) == prob

    frac = 0.1
    rule = WaterUse(spp, loc, frac)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
    @test rule.soil_moisture_frac == frac

    rule = WindDispersal(spp, loc, 1.0m, 0.3, 1.0m, 1.0m/s, 1.0m/s)
    @test getspecies(rule) == spp
    @test getlocation(rule) == loc
end