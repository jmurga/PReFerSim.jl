################################
####### Define parameters ######
################################
@with_kw struct recipe
    epochs::Array{Int64,1} = [1000]
    N::Array{Int64,1} = [10000]

    θ::Float64=8
    h::Float64=0.5
    
    s₋::Float64=-457.0
    s₊::Float64=500.0

    dfe::String="point"
    param_one::Float64=1.0
    param_two::Float64=1.0
    s::Array{Float64,1}=-[0.0]
    s_mult::Array{Float64,1}=[1.0]
    prob::Array{Float64,1}=[0.0]

    n_anc::Int64=1000
    burnin_period::Bool=false
    
    relax::Bool=false
    epoch_relaxation::Array{Bool}=fill(false,length(epochs))
    s_relaxation::Float64=0.0
    s_relaxation_threshold::Float64=0.0

    F::Array{Float64,1}=zeros(length(N))
    seed::Int64=rand(1:10^8)

    trajectories::Array{Int64,1} = Int64[]
    trajectories_output::String = ""

    @assert length(N)==length(epochs)  "N and epochs must be equal in length";
    @assert length(s)==length(prob)    "s and probs must be equal in length";
end

@with_kw mutable struct prf_output
    fixations::Int64 = 0;
    loss::Int64=0;
    total_mut::Int64 = 0;
    count_mut::Int64=0;
end

@with_kw mutable struct mutation
    # The frequency of the mutation
    frequency::Float64=0;
    # Count of the mutation
    count_samp::Int64=1; 
    # Selection coefficient
    s::Float64=0.0; 
    # Dominance factor
    h::Float64=0.5; 
    # The generation that mutation has arise
    age::Int64=1;
    # Number of mutation. Tag to avoid similar nodes
    num::Int64=1;
    # Mutation type. 1=neutral; 2=deleterious; 3=strong advantegeous; 4=weakly advantegeous
    type::Int8=0
end
