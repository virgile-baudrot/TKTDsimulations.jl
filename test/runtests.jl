using TKTDsimulations
using Test
using DifferentialEquations

@testset "Linear interpolation" begin
    # TEST INTERNAL FUNCTIONS: need specified TKTDsimulations because not Export
    # 1. linear interpolation
    @test TKTDsimulations.interpLinear([0,1,2,3], [0,0,2,2], 1.5) == 1.0 
    @test TKTDsimulations.interpLinear([0,1,2,3], [0,0,2,2], 5.0) == 0
end

@testset "TK solver" begin
    # 2. odeTK
    @test_nowarn begin
        probTK = ODEProblem(TKTDsimulations.odeTK, [0.0], (0.0,3.0), [[0,1,2,3], [0,0,2,2], 0.5])
        _saveatTK = [0,1,2,3] === nothing ? Float64[] : [0,1,2,3]
        solve(probTK ; saveat = _saveatTK)
    end

    # tps = [0,1,2,3]
    # conc = [0,0,2,2]
    # kd = 0.5
    @test_nowarn runTK([0,1,2,3], [0,0,2,2], 0.5)
    # kd = [0.5, 0.6]
    @test_nowarn runTK_MCMC([0,1,2,3], [0,0,2,2], [0.5, 0.6])
end

@testset "SD solver" begin
    # 2.2 odeSD
    @test_nowarn begin
        probSD = ODEProblem(TKTDsimulations.odeSD, [0.0, 0.0], (0.0,3.0), [[0,1,2,3], [0,0,2,2], 0.5, 0.2, 0.3, 0.4])
        _saveatSD = [0,1,2,3] === nothing ? Float64[] : [0,1,2,3]
        solve(probSD ; saveat = _saveatSD)
    end

    tps = [0,1,2,3]
    conc = [5,5,5,5]
    kd = 0.5
    hb = 0.5
    z = 0.5
    kk = 0.5
    # TEST EXPORT FUNCTION: no need of TKTDsimulations specification because Export
    @test_nowarn runSD(tps, conc, kd, hb, z, kk)

    kd = [0.5, 2.0]
    hb = [0.5, 2.0]
    z = [0.5, 2.0]
    kk = [0.5, 2.0]
    @test_nowarn runSD_MCMC(tps, conc, kd, hb, z, kk)

end

@testset "IT solver" begin

    tps = [0,1,2,3]
    conc = [0,0,2,2]
    kd = 0.5
    hb = 0.5
    alpha = 0.5
    beta = 0.5
    @test_nowarn runIT(tps, conc, kd, hb, alpha, beta)

    kd = [0.5, 2.0]
    hb = [0.5, 2.0]
    alpha = [0.5, 2.0]
    beta = [0.5, 2.0]
    @test_nowarn runIT_MCMC(tps, conc, kd, hb, alpha, beta)

end
