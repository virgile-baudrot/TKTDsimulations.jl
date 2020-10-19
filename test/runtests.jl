using TKTDsimulations
using Test
using DifferentialEquations

@testset "TKTDsimulations.jl" begin
    # TEST INTERNAL FUNCTIONS: need specified TKTDsimulations because not Export
    # 1. linear interpolation
    @test TKTDsimulations.interpLinear([0,1,2,3], [0,0,2,2], 1.5) == 1.0 
    @test TKTDsimulations.interpLinear([0,1,2,3], [0,0,2,2], 5.0) == 0

    # 2. ode solver 
    # 2.1 odeTK
    @test_nowarn begin
        probTK = ODEProblem(TKTDsimulations.odeTK, [0.0], (0.0,3.0), [[0,1,2,3], [0,0,2,2], 0.5])
        _saveatTK = [0,1,2,3] === nothing ? Float64[] : [0,1,2,3]
        solve(probTK ; saveat = _saveatTK)
    end
    # 2.2 odeSD
    @test_nowarn begin
        probSD = ODEProblem(TKTDsimulations.odeSD, [0.0, 0.0], (0.0,3.0), [[0,1,2,3], [0,0,2,2], 0.5, 0.2, 0.3, 0.4])
        _saveatSD = [0,1,2,3] === nothing ? Float64[] : [0,1,2,3]
        solve(probSD ; saveat = _saveatSD)
    end

    # TEST EXPORT FUNCTION: no need of TKTDsimulations specification because Export
    @test_nowarn runTK([0,1,2,3], [0,0,2,2], 0.5)
    @test_nowarn runSD([0,1,2,3], [0,0,2,2], 0.5, 0.2, 0.3, 0.4)
    @test_nowarn runIT([0,1,2,3], [0,0,2,2], 0.5, 0.2, 0.3, 0.4)

end
