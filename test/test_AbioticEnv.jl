using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using Peatland
using AxisArrays
import EcoSISTEM: getlocation

grid = (5, 5)
area = 25.0km^2
totalK = 10000.0kJ/km^2
numNiches = 4
active = fill(true, grid)

@testset "peatland habitat" begin
    # TEST peatlandAE
    fillval = 10.0m^3
    abenv = peatlandAE(fillval, grid, totalK, area)
    @test_nowarn peatlandAE(fillval, grid, totalK, area)
    @test_nowarn peatlandAE(fillval, grid, totalK, area, active)
    @test all(abenv.habitat.matrix .== fillval)
    @test size(abenv.habitat.matrix) == grid
    @test sum(abenv.budget.matrix) == totalK * area
    @test abenv.active == active
    @test all(abenv.active)

    @test EcoSISTEM._getsubcommunitynames(abenv) == abenv.names
    @test EcoSISTEM.getavailableenergy(abenv) == sum(abenv.budget.matrix)
end

@testset "elevation peatland habitat" begin
    # TEST peatlandAE
    fillval = 10.0m^3
    ele = AxisArray(fill(1.0, grid), Axis{:latitude}(1.0km:1.0km:5km), Axis{:longitude}(1.0km:1.0km:5km))
    abenv = peatlandAE(ele, fillval, totalK)
    @test_nowarn peatlandAE(ele, fillval, totalK)
    @test_nowarn peatlandAE(ele, fillval, totalK, active)
    @test all(abenv.habitat.h1.matrix .== fillval)
    @test size(abenv.habitat.h1.matrix) == grid
    @test all(abenv.habitat.h2.matrix .== 1.0)
    @test size(abenv.habitat.h2.matrix) == grid
    @test sum(abenv.budget.matrix) == totalK * area
    @test abenv.active == active
    @test all(abenv.active)

    @test EcoSISTEM._getsubcommunitynames(abenv) == abenv.names
    @test EcoSISTEM.getavailableenergy(abenv) == sum(abenv.budget.matrix)

    loc = 1; prob = 0.1/day
    rule = LateralFlow(abenv, loc, prob)
    @test getlocation(rule) == loc
    rule = LateralFlow(abenv, abenv.habitat.h2, loc, prob)
    @test getlocation(rule) == loc
end

