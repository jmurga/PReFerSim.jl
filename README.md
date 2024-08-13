# PReFerSim

[![Build Status](https://github.com/jmurga/PReFerSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jmurga/PReFerSim.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a Julia port of PReFerSim [PReFerSim](https://doi.org/10.1093/bioinformatics/btw478) software. It simplify parameters input, modification and simulation replicas. PReFerSim performs forward simulations under the PRF model while allowing models changes in population size, arbitrary amounts of inbreeding, dominance, and distributions of selective effects.

If you use PReFerSim please cite: Ortega-Del Vecchyo, D., Marsden, C. D. & Lohmueller, K. E. PReFerSim: fast simulation of demography and selection under the Poisson Random Field model. Bioinformatics 32, 3516–3518 (2016).

The package is follow same parameter inputs as the orginal PReFerSim implementation. Check the original  manual for a comprehensive overview of the input parameters.

Now the parameteres are accesible by using `PreFerSim.recipe()` structure. The following example mimics ParameterFile1.txt

```julia
using PReFerSim

demes = PReFerSim.recipe(N = [18169,44993,1999,2000,1000],epochs = [145352,63898,347,2033,100], θ = 8)
sample_size = 20
sfs, fix = PReFerSim.simulate(demes,sample_size)
```

Multiple recipes can be runned at once by taking advantage of of Julia multithreading (you must start Julia by using `-t<nthreads>`;e.g `julia -t8`)


```julia
using PReFerSim

demes = PReFerSim.recipe(N = [18169,44993,1999,2000,1000],epochs = [145352,63898,347,2033,100], θ = 8)

# Custom Vector containing 10 times the previous recipe. You can pool the results by using PReFerSim.simulate(pool=true)

demes_replicas = fill(demes,10);
sample_size = 20
sfs, fix = PReFerSim.simulate(demes_replicas,sample_size)
```
