function fs(mutation_list::Vector{Mutation}, sample_size::Int, r::Ptr{gsl_rng})
    # Pre-allocate output: column 1 = bin center/label, column 2 = counts
    out = zeros(sample_size - 1, 2)

    @inbounds for i in 1:(sample_size - 1)
        out[i, 1] = round(i / sample_size, digits = 3)
    end

    isempty(mutation_list) && return out

    invS = 1.0 / sample_size

    @inbounds for m in mutation_list
        freq = m.frequency

        # GSL binomial draw; trials = sample_size
        samp_count = ran_binomial(r, freq, sample_size)
        samp_freq  = samp_count * invS

        # Count only polymorphic samples
        if 0.0 < samp_freq < 1.0
            # Map samp_count ∈ {1, …, sample_size-1} to bin index ∈ {1, …, sample_size-1}
            bin_idx = samp_count

            out[bin_idx, 2] += 1.0
        end
    end

    return out
end


@inline loginfo(::Val{true},  msg) = (@info msg; nothing)
@inline loginfo(::Val{false}, msg) = nothing


@inline logwarn(::Val{true},  msg) = (@warn msg; nothing)
@inline logwarn(::Val{false}, msg) = nothing


function simulate(param::recipe, sample_size::Int64;verbose=true)
    
    V = Val(verbose)

    @unpack epochs, N,θ,h,s₋,s₊,dfe,param_one,param_two,s,s_mult,prob,n_anc,burnin_period,relax,epoch_relaxation,s_relaxation,s_relaxation_threshold,F,trajectories,seed = param
    
    # Create RNG for this simulation (thread-safe)
    T = gsl_rng_default
    r = rng_alloc(T)
    rng_set(r, seed)
    
    # Create state for trajectory tracking
    state = isempty(trajectories) ? nothing : Trajectories(Int32(0),Dict(Int32(id) => Float64[] for id in trajectories))

    # Set up DFE distribution
    if dfe == "point"
        dfe_dist = s
    elseif dfe == "gamma"
        # Use first s_mult value (or could be epoch-specific later)
        dfe_dist = Gamma(param_one, param_two * s_mult[1])
    elseif dfe == "beta"
        dfe_dist = Beta(param_one, param_two * s_mult[1])
    elseif dfe == "lognormal"
        dfe_dist = LogNormal(param_one, param_two * s_mult[1])
    end

    events::Int64 = length(epochs)
    
    loginfo(V,"Demographic History ($(events) epochs)")
    for e = 1:events
        #@info "Ne = $(N[e]) generations = $(epochs[e]); F = $(F[e])"
        loginfo(V, "Ne = $(N[e]) generations = $(epochs[e]); F = $(F[e])")
    end
    
    # Burnin
    if burnin_period && dfe != "point"
        #@warn "Burnin period will assume neutrality  (dfe == point; s = 0"
        logwarn(V, "Burnin period will assume neutrality (dfe == point; s = 0)")
    end

    mutation_list = burnin_period ? burnin(param, [0], r) : Vector{Mutation}()

    l::Int64 = 0
    f::Int64 = 0
    age::Int64 = 0
    θ_dist = Poisson(θ / 2.0)


    @inbounds for e = 1:events
        # Inbreeding Ne

        N_F::Float64 = N[e] / (1.0 + F[e]);

        # @info "Currently in epoch = $(e) ; Segregating mutations = $(length(mutation_list))"
        loginfo(V, "Currently in epoch = $e ; Segregating mutations = $(length(mutation_list))")

        if epoch_relaxation[e]
            # @info "Relaxation in Epoch $(e)"
            loginfo(V, "Relaxation in Epoch $e")
            relax_selection!(mutation_list, s_relaxation, s_relaxation_threshold, 0)
        end
        
        sizehint!(mutation_list, length(mutation_list) + Int(ceil(mean(θ_dist) * epochs[e])))

        N_freq = 1.0 / N_F
        for g::Int64 = 1:epochs[e]
            # Mutation age
            age += 1
            
            # Drift and selection
            x, y = drift_sel!(mutation_list, r, N_F, F[e], state)
            
            # Add new mutations
            add_mutation!(
                mutation_list,
                N_F,
                h,
                θ_dist,
                N_freq,
                dfe_dist,
                n_anc,
                age,
                state,
                r
            )
            
            l += x
            f += y
        end
        

        # Update theta for next epoch
        if e < events
            N_F_1 = N[e + 1] / (1.0 + F[e + 1])
            θ = θ * N_F_1 / N_F
            θ_dist = Poisson(θ / 2.0)
        end
    end
    
    # Compute SFS using same RNG
    out = fs(mutation_list, sample_size, r)
    
    #@info "Total number of mutations = $(l + f)"
    loginfo(V, "Total number of mutations = $(l + f)")

    
    # Return results
    if state === nothing
        return (out, hcat(l, f))
    else
        return (out, hcat(l, f), state.trajectories_output)
    end

end


function simulate(param::Vector{recipe},sample_size::Int64;pool::Bool=false)

    @info "Running a total of $(length(param)) recipes in $(Threads.nthreads()) threads"
    
    sfs,fix = @suppress begin
        unzip(ThreadsX.map(x -> simulate(x,sample_size,verbose=false),param));
    end

    if pool
        sfs = sum(sfs);
        sfs[:,1] .= round.(collect(1:(sample_size-1))/sample_size,digits=3)
        fix = sum(fix)
    end

    return sfs,fix
end

