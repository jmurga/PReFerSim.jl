function fs(mutation_list::LinkedList{mutation},r::Ptr{gsl_rng},sample_size::Int64)

    if(!isempty(mutation_list))

        mut  = zeros(length(mutation_list)); 
        #=neut = zeros(length(mutation_list));
        del  = zeros(length(mutation_list));
        adv  = zeros(length(mutation_list));=#
        
        i    = 1;
        # Index to delete node
        node = mutation_list.node.next;

        while (node.next != mutation_list.node.next)
            item        = node.data;

            freq        = item.frequency;
            s           = item.s;

            h           = item.h;

            #samp_count = rand(Binomial(sample_size,freq));
            samp_count = ran_binomial(r,freq,sample_size);
            samp_freq = samp_count/sample_size;
            if (samp_freq > 0 && samp_freq < 1.0)
                mut[i] = samp_freq
                #=if(item.type == 1)
                    neut[i] = samp_freq
                elseif(item.type == 2)
                    del[i] = samp_freq
                elseif(item.type == 3)
                    adv[i] = samp_freq
                end=#
            end
            node = node.next;
            i+=1
        end

        mut     = mut[mut .!= 0]

        b        = collect(1/sample_size:1/sample_size:1)
        v(x,b=b) = searchsortedfirst.(Ref(b), x)
        n    = v(mut)
        
        out      = zeros((sample_size-1,2))
        out[:,1] = collect(1:(sample_size-1))

        @inbounds for i::Int64 ∈ out[:,1]
            out[i,2] = length(n[n .== i])
        end

        out[:,1] = round.(out[:,1]/sample_size,digits=3);

        return out
    else
        out = zeros((sample_size-1,3));
        out[:,1] = round.(out[:,1]/sample_size,digits=3);
        return out
    end
end

function simulate(param::recipe, sample_size::Int64)
    
    @unpack epochs,N,θ,h,s₋,s₊,dfe,param_one,param_two,s,s_mult,prob,n_anc,burnin_period,relax,epoch_relaxation,s_relaxation,s_relaxation_threshold,F,trajectories,trajectories_output,seed = param;
    
    seed = rand(1:10^8)
    if !isempty(trajectories)
        @assert length(unique(trajectories .> 0)) == 1  "ID index must be greater than 0";
        # global io = open(trajectories_output,"w+")

        if isempty(trajectories_output)
            trajectories_output = OrderedDict(trajectories .=> Vector{Float64}[[]])
        end

        global trajectories_output
    end
    
    # set seed before start
    T = gsl_rng_default;
    r = rng_alloc(T);

    rng_set(r, seed)

    epochs = SVector{length(epochs)}(epochs)
    s_mult = SVector{length(s_mult)}(s_mult)
    s = SVector{length(s)}(s)
    F = SVector{length(F)}(F)

    events::Int64 = length(epochs);
    
    if(length(s)==1)
        s      = s[1];
        s_mult = s_mult[1]
    else
        s = s;
        s_mult = s_mult;
    end
    

    @info "Demographic History ($(events) epochs)"
    for e=1:events
        @info "Ne = $(N[e]) generations = $(epochs[e]) F = $(param.F[e])"
    end

    if burnin_period
        mutation_list = burnin(param,r);
    else
        mutation_list = LinkedList{mutation}();
    end;
    
    l::Int64= 0;
    f::Int64= 0;
    N_F_1::Int64=0;
    # @printf "\n"
    @inbounds for e=1:events
        # Inbreeding Ne

        N_F::Float64 = N[e] / (1.0 + F[e]);

        @info "Currently in epoch = $(e) ; Segregating mutations before epoch's beginning = $(length(mutation_list))";


        if epoch_relaxation[e]
            @info "Relaxation in Epoch $(e)";

            relax_selection!(mutation_list, s_relaxation, s_relaxation_threshold, relaxation_type);
        end

        @inbounds for g::Int64=1:epochs[e]
            # Mutation age
            if e == 1 
                age = g;
            else
                # Add age to previous epochs generations
                age = g + epochs[e];
            end
            x, y    = drift_sel!(mutation_list,r,N_F,F[e],trajectories);
            add_mutation!(mutation_list,r,N_F,h,s,θ,1.0/N[e],dfe,param_one,param_two,s_mult,n_anc,age,trajectories,relax);
            l = l + x;
            f = f + y;
        end

        if (e < events) 
            N_F_1 = N[e + 1] / (1.0 + param.F[e + 1]);
            θ = θ * N_F_1 / N_F;
        end
    end

    # if @isdefined io
    #     close(io)
    # end

    out = fs(mutation_list,r,sample_size);

    @info "Total number of mutations = $(l+f)"


    if isempty(trajectories_output)
        return(out,hcat(l,f))
    else
        return(out,hcat(l,f),trajectories_output)
    end
end

function simulate(param::Vector{recipe},sample_size::Int64;pool::Bool=false)

    @info "Running a total of $(length(param)) recipes in $(Threads.nthreads())"
    sfs,fix = @suppress begin
        unzip(ThreadsX.map(x -> simulate(x,sample_size),param));
    end

    if pool
        sfs = sum(sfs);
        sfs[:,1] .= round.(collect(1:(sample_size-1))/sample_size,digits=3)
        fix = sum(fix)
    end
    return sfs,fix
end

