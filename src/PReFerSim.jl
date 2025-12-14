module PReFerSim

using GSL, Parameters, QuadGK, ThreadsX

import Distributions: Gamma, Beta, Poisson, Binomial, LogNormal, pdf, mean
import OrderedCollections: OrderedDict
import Unzip:unzip

include("parameters.jl")
include("prf_drift.jl")
include("simulate.jl")

end
