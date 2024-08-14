using PReFerSim
using Test

@testset "PReFerSim.jl" begin

    demes = PReFerSim.recipe(trajectories = [1,5,10],burnin_period=true)
    @test demes.epochs[1] == 200000

end
