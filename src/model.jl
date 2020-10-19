function interpLinear(xpt,ypt,x)
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
    du[2] = p[4] + maximum([u[2] - p[5], 0]) * p[6]
end

function logLogisticLaw(x, α, β)
    1 /(1 +(x/α)^(-β))
end

"""
Solve ODE for Stochastic Death TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `p` -- parameters
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
Solve ODE for Stochastic Death TKTD model

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- parameters
- `z` -- parameters
- `kk` -- parameters
"""
function runSD(tps, conc, kd, hb, z, kk)
    prob = DiffEqBase.ODEProblem(odeSD,[0.0, 0.0],(0.0, maximum(tps)),[tps, conc, kd, hb, z, kk])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)
    return DataFrames.DataFrame(
        time = tps,
        exposure = conc,
        TK = [sol.u[j][1] for j in 1:length(sol)],
        TD = [sol.u[j][2] for j in 1:length(sol)]
    )
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
    for j in 1:length(Dmax_tmp)
        pSurv[j] = exp(-hb*tps[j])*(1-logLogisticLaw(Dmax_tmp[j][1], alpha,beta))
    end
    return DataFrames.DataFrame(
        time = tps,
        exposure = conc,
        TK = [sol.u[j][1] for j in 1:length(sol)],
        TD = pSurv
    )
end