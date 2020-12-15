# TKTDsimulations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://virgile-baudrot.github.io/TKTDsimulations.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://virgile-baudrot.github.io/TKTDsimulations.jl/dev)
[![Build Status](https://github.com/virgile-baudrot/TKTDsimulations.jl/workflows/CI/badge.svg)](https://github.com/virgile-baudrot/TKTDsimulations.jl/actions)
[![Coverage](https://codecov.io/gh/virgile-baudrot/TKTDsimulations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/virgile-baudrot/TKTDsimulations.jl)

## Examples

```julia
julia> using TKTDsimulations

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

#### For registration

Registration is done using [Registrator.jl](https://github.com/JuliaRegistries/Registrator.jl).

For the registration of this package, I did it with *Register pacakges* in [https://juliahub.com](https://juliahub.com)

For its updates: I open the issue [Registrator trigger](https://github.com/virgile-baudrot/TKTDsimulations.jl/issues/8) and I add comment mentioning `@JuliaRegistrator register` and a realease notes (maybe more elaborate compared to these first release)

```julia
@JuliaRegistrator register

Release notes:

Check out my new features!
```

Then,  a `Pull Request` is done from the Julia Tagbot, and once merged the package is updated in JuliaHub. See the issues comments [Registrator trigger](https://github.com/virgile-baudrot/TKTDsimulations.jl/issues/8) for information. 

After the above pull request is merged, it is recommended that a tag is created on this repository for the registered package version. Or, this will be done automatically if the Julia TagBot GitHub Action is installed.
