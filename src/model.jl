function interpLinear(xpt::AbstractArray,ypt::AbstractArray,x)
    if x >= minimum(xpt) && x <= maximum(xpt)
        idx = findfirst(t->t>= x, xpt) -1
        idx == 0 ? idx += 1 : nothing
        θ = (x - xpt[idx])/ (xpt[idx+1] - xpt[idx])
        (1-θ)*ypt[idx] + θ*ypt[idx+1]
    else
        0
    end
end

function odeTK(du, u, p,t)
    du[1] = p[3] * (interpLinear(p[1],p[2],t)-u[1])
end

function odeSD(du, u, p,t)
    du[1] = p[3] * (interpLinear(p[1],p[2],t)-u[1])
    du[2] = p[4] + maximum([u[1] - p[5], 0]) * p[6]
end

function logLogisticLaw(x, α, β)
    1 /(1 +(x/α)^(-β))
end

"""
    Solve Toxicokinetics Models

**Fields**

    - `tps` -- time vector
    - `conc` -- exposure vector
    - `kd` -- parameter, scalar
"""
function runTK(tps, conc, kd)
    prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)
    return DataFrames.DataFrame(
        time = tps,
        exposure = conc,
        TK = [sol.u[j][1] for j in 1:length(sol)]
    )
end

"""
    Solve Toxicokinetics Models with MCMC parameters

**Fields**

    - `tps` -- time vector
    - `conc` -- exposure vector
    - `kd` -- parameter, scalar
"""
function runTK_MCMC(tps, conc, kd)
    TK = Array{Array{Array{Float64,1},1}, 1}(undef, length(kd))
    @inbounds for i in eachindex(kd)
        prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd[i]])
        _saveat = tps === nothing ? Float64[] : tps
        TK[i] = DiffEqBase.solve(prob ; saveat = _saveat).u
    end
    return TK 
end


"""
Solve ODE for Stochastic Death TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- background mortality
- `z` -- threshold
- `kk` -- killing rate
"""
function runSD(tps, conc, kd, hb, z, kk)
    prob = DiffEqBase.ODEProblem(odeSD,[0.0, 0.0],(0.0, maximum(tps)),[tps, conc, kd, hb, z, kk])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)
    TK = [sol.u[j][1] for j in 1:length(sol)]
    TD = [exp(-sol.u[j][2]) for j in 1:length(sol)]
    return TK, TD
end


"""
Solve ODE for Stochastic Death TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- background mortality
- `z` -- threshold
- `kk` -- killing rate
"""

function runSD_MCMC(tps, conc, kd, hb, z, kk)
    TK = Array{Array{Float64,1},1}(undef, length(kd))
    TD = Array{Array{Float64,1},1}(undef, length(kd))
    @inbounds for i in eachindex(kd)
        prob = DiffEqBase.ODEProblem(odeSD,[0.0, 0.0],(0.0, maximum(tps)),[tps, conc, kd[i], hb[i], z[i], kk[i]])
        _saveat = tps === nothing ? Float64[] : tps
        sol = DiffEqBase.solve(prob ; saveat = _saveat)
        TK[i] = [sol.u[j][1] for j in 1:length(sol)]
        TD[i] = [exp(- sol.u[j][2]) for j in 1:length(sol)]
    end
    return TK, TD
end


"""
Solve ODE for Individual Tolerance TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- parameters
- `alpha` -- parameters
- `beta` -- parameters
"""
function runIT(tps, conc, kd, hb, alpha, beta)
    prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)

    pSurv = Array{Float64}(undef,length(tps))
    Dmax_tmp = accumulate(max, sol.u)

    TK = [sol.u[j][1] for j in 1:length(sol)]
    Dmax_tmp = accumulate(max, sol.u)      
    TD = [ exp(-hb*tps[j])*(1-logLogisticLaw(Dmax_tmp[j][1], alpha,beta)) for j in 1:length(sol)]
    
    return TK, TD
end


"""
Solve ODE for Individual Tolerance TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- parameters
- `alpha` -- parameters
- `beta` -- parameters
"""

function runIT_MCMC(tps, conc, kd, hb, alpha, beta)
    TK = Array{Array{Float64,1},1}(undef, length(kd))
    TD = Array{Array{Float64,1},1}(undef, length(kd))
    @inbounds for i in eachindex(kd)
        prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd[i]])
        _saveat = tps === nothing ? Float64[] : tps
        sol = DiffEqBase.solve(prob ; saveat = _saveat)
        TK[i] = [sol.u[j][1] for j in 1:length(sol)]
        Dmax_tmp = accumulate(max, sol.u)      
        TD[i] = [ exp(-hb[i]*tps[j])*(1-logLogisticLaw(Dmax_tmp[j][1], alpha[i],beta[i])) for j in 1:length(sol)]
    end
    return TK, TD
end


# struct ExposureProfile{T<:Number,X<:Number}
#     time::Array{T,1}
#     exposure::Array{X,1}
#     function ExposureProfile(time::Array{T,1},
#                              exposure::Array{X,1}) where {T,X}
#         if length(time) != length(exposure)
#             throw(ArgumentError("length of time and exposure differ"))
#         end
#         if length(time) != length(unique(time))
#             throw(ArgumentError("values in time vector must be unique"))
#         end
#         if sort(time) != time
#             throw(ArgumentError("time vector must be sorted"))
#         end
#         new{T,X}(time, exposure)
#     end
# end
# ExposureProfile([1.0],[2.0])

# ExposureProfile([1.0,2,3],[2.0,3.0,2.0])
# ExposureProfile([1,2,3],[2.0,3.0,2.0])
# ExposureProfile([1,2,3],[2.0,3.0,2.0]).time
# ExposureProfile([1,2,3],[2.0,3.0,2.0]).exposure
# # should return false
# ExposureProfile([1,3,3],[2.0,3.0,2.0])
# ExposureProfile([1,3,2],[2.0,3.0,2.0])
# ExposureProfile([1.0],[2.0,3.0])