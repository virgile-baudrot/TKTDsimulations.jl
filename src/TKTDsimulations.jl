module TKTDsimulations

    import DiffEqBase
    import DifferentialEquations
    import DataFrames

    include("model.jl") 

    export runTK
    export runSD
    export runIT
    # export runGUTS

end
