module PReFerSim

using GSL, LinkedLists, Parameters, StaticArrays, Printf, QuadGK, ProgressMeter, ThreadsX

include("parameters.jl")
include("prf_drift.jl")
include("simulate.jl")

end
