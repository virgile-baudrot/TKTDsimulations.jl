module TKTDsimulations

    import DiffEqBase
    import DifferentialEquations

    include("model.jl") 

    export runTK, runTK_MCMC
    export runSD, runSD_MCMC
    export runIT, runIT_MCMC

end
