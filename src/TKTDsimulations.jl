module TKTDsimulations

    import DiffEqBase
    import DifferentialEquations
    # import DataFrames

    include("model.jl") 

    export runTK, runTK_MCMC
    export runSD, runSD_MCMC
    export runIT, runIT_MCMC
    # export runGUTS

end
