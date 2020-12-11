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
    Solve Toxicokinetics Models with a single parameter

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

**Fields**
    - `tps` -- time vector
    - `conc` -- exposure vector
    - `kd` -- parameter, a scalar

Return a Tuple with `Array{Float64,1}` of th same length as `tps`

# Example
```julia-repl
julia> myTK = runTK([0,1,2,3], [0,1,20,2], 0.5);
julia> myTK.TK
4-element Array{Float64,1}:
 0.0
 0.21800446293645104
 4.570824232907027
 6.801199145490319
```
"""
function runTK(tps, conc, kd)
    prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)
    TK = [sol.u[j][1] for j in 1:length(sol)]
    return (TK=TK,)
end

"""
    Solve Toxicokinetics Models with MCMC parameters

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

**Fields**

    - `tps` -- time vector
    - `conc` -- exposure vector
    - `kd` -- parameter, scalar

Return an `Array{Array{Float64,1},1}`.

# Example
```julia-repl
julia> myTK = runTK_MCMC([0,1,2,3], [0,1,20,2], [0.5, 2.4]);

julia> myTK.TK
2-element Array{Array{Float64,1},1}:
 [0.0, 0.21800446293645104, 4.570824232907027, 6.801199145490319]
 [0.0, 0.6019164649918195, 12.77507909562312, 8.171280475528322]
```
"""
function runTK_MCMC(tps, conc, kd)
    TK = Array{Array{Float64,1},1}(undef, length(kd))
    @inbounds for i in eachindex(kd)
        prob = DiffEqBase.ODEProblem(odeTK,[0.0],(0.0, maximum(tps)),[tps, conc, kd[i]])
        _saveat = tps === nothing ? Float64[] : tps
        sol = DiffEqBase.solve(prob ; saveat = _saveat)
        TK[i] = [sol.u[j][1] for j in 1:length(sol)]
    end
    return (TK=TK,)
end


"""
Solve ODE for Stochastic Death TKTD model

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

    \$\\frac{dy(t)}{dt} = k_d (x(t) - z) + h_b\$

    \$ pSurv(t) = \\exp{(-y(t))}\$

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- background mortality
- `z` -- threshold
- `kk` -- killing rate

Return a Tuple of `Array{Float64,1}`.

# Example
```julia-repl
julia> mySD = runSD([0,1,2,3], [0,1,20,2], 0.5, 2.4, 16.3 ,10.0);

julia> mySD.TK
4-element Array{Float64,1}:
 0.0
 0.21089023114020433
 4.575797874578747
 6.801834882676336

julia> mySD.TD
4-element Array{Float64,1}:
 1.0
 0.09071795328941247
 0.008229747049020037
 0.0007465858083766806
```
"""
function runSD(tps, conc, kd, hb, z, kk)
    prob = DiffEqBase.ODEProblem(odeSD,[0.0, 0.0],(0.0, maximum(tps)),[tps, conc, kd, hb, z, kk])
    _saveat = tps === nothing ? Float64[] : tps
    sol = DiffEqBase.solve(prob ; saveat = _saveat)
    TK = [sol.u[j][1] for j in 1:length(sol)]
    TD = [exp(-sol.u[j][2]) for j in 1:length(sol)]
    return (TK=TK, TD=TD,)
end


"""
Solve ODE for Stochastic Death TKTD model

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

    \$\\frac{dy(t)}{dt} = k_d (x(t) - z) + h_b\$

    \$ pSurv(t) = \\exp{(-y(t))}\$

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- parameters
- `hb` -- background mortality
- `z` -- threshold
- `kk` -- killing rate

Return a Tuple of `Array{Array{Float64,1},1}`.

# Example
```julia-repl
julia> mySD = runSD_MCMC([0,1,2,3], [0,1,20,2], [0.5, 2.4],[2.5, 2.1],[10.5, 16.0],[5.0, 2.4]);

julia> mySD.TK
2-element Array{Array{Float64,1},1}:
 [0.0, 0.21089023114020433, 4.575797874578747, 6.801834882676336]
 [0.0, 0.6309497946830392, 12.689036317855038, 8.153727962011798]

julia> mySD.TD
2-element Array{Array{Float64,1},1}:
 [1.0, 0.08208499862389876, 0.006737946999085473, 0.0005530843701478346]
 [1.0, 0.12245642825298195, 0.014995576820477731, 0.0018363047770289104]
```

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
    return (TK=TK, TD=TD)
end


"""
Solve ODE for Individual Tolerance TKTD model

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

    \$ pSurv(t) = \\exp{(-h_b t)} (1-logLogistic(max(x(t), alpha,beta))  \$

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- toxicokinetics parameter
- `hb` -- background mortality
- `alpha` -- median of log-logistic distribution
- `beta` -- shape of the log-logistic distribution

Return an `Array{Float64,1}`.

# Example
```julia-repl
julia> myIT = runIT([0,1,2,3], [0,1,20,2], 0.5, 2.1, 10.5 , 2.4);

julia> myIT.TK
4-element Array{Float64,1}:
 0.0
 0.21800446293645104
 4.570824232907027
 6.801199145490319
 
julia> myIT.TD
4-element Array{Float64,1}:
 1.0
 0.12244522344410577
 0.013201809910381767
 0.001357555235565695
```

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

    return (TK=TK, TD=TD)
end


"""
Solve ODE for Individual Tolerance TKTD model

    \$\\frac{dx(t)}{dt} = k_d (conc(t) - x(t))\$

    \$ pSurv(t) = \\exp{(-h_b t)} (1-logLogistic(max(x(t), alpha,beta))  \$

**Fields**

- `tps` -- time vector
- `conc` -- exposure vector
- `kd` -- toxicokinetics parameter
- `hb` -- background mortality
- `alpha` -- median of log-logistic distribution
- `beta` -- shape of the log-logistic distribution

Return an `Array{Array{Float64,1},1}`.

# Example
```julia-repl
julia> myIT = runIT_MCMC([0,1,2,3], [0,1,20,2], [0.5, 2.4],[2.5, 2.1],[10.5, 16.0],[5.0, 2.4]);

julia> myIT.TK
2-element Array{Array{Float64,1},1}:
 [0.0, 0.21800446293645104, 4.570824232907027, 6.801199145490319]
 [0.0, 0.6019164649918195, 12.77507909562312, 8.171280475528322]

julia> myIT.TD
2-element Array{Array{Float64,1},1}:
 [1.0, 0.0820849983072016, 0.006634237831805787, 0.0004964761672208183]
 [1.0, 0.12240978217836299, 0.009475165156808458, 0.0011602948822098684]
```

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
    return (TK=TK, TD=TD)
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