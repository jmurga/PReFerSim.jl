################################
####### Define parameters ######
################################
@with_kw mutable struct recipe

    epochs::Vector{Int64} = [200000]
    N::Vector{Int64} = [10000]

    θ::Float64=8
    h::Float64=0.5
    
    s₋::Float64=-457.0
    s₊::Float64=500.0

    dfe::String="point"
    param_one::Float64=1.0
    param_two::Float64=1.0

    s::Vector{Float64}=-[0.0]
    s_mult::Vector{Float64}=[1.0]
    prob::Vector{Float64}=[0.0]

    # n_anc = param_three, used as oldest N
    n_anc::Int64=N[1]
    burnin_period::Bool=false
    
    relax::Bool=false
    epoch_relaxation::Array{Bool}=fill(false,length(epochs))
    s_relaxation::Float64=0.0
    s_relaxation_threshold::Float64=0.0

    F::Vector{Float64}=zeros(length(N))

    trajectories::Vector{Int64} = Int64[]

    seed::Int64 = rand(1:10^8)
    #trajectories_output::OrderedDict{Int64,Vector} = ifelse(isempty(trajectories),OrderedDict{Int64,Vector}(),OrderedDict{Int64,Vector}(trajectories .=> Vector{Float64}[[]]))

    @assert length(N)==length(epochs)  "N and epochs must be equal in length";
    @assert length(s)==length(prob)    "s and probs must be equal in length";
end

@with_kw mutable struct prf_output
    fixations::Int64 = 0;
    loss::Int64=0;
    total_mut::Int64 = 0;
    count_mut::Int64=0;
end

struct Mutation
    frequency::Float64
    h::Float32
    s::Float64
    count_samp::UInt32
    age::UInt32
    num::UInt32
    type::UInt8
end

mutable struct Trajectories
    next_mut_id::Int32
    trajectories_output::Dict{Int32, Vector{Float64}}
end