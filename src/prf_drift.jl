function relax_selection!(mutation_list::LinkedList{mutation}, new_s, sel_threshold::Float64, relax_type::Int64)

    if(!isempty(mutation_list))

        # Index to delete node
        node=mutation_list.node.next;

        while (node.next != mutation_list.node.next)
            item        = node.data;

            freq            = item.frequency;
            s               = item.s;
            h               = item.h;
            freq            = freq_inbreed(s,freq,h);
            num       = item.num;

            if ((s/2) <= sel_threshold && relax_tye == 0)
                item.s = new_s * 2;
            elseif ((s/2) <= sel_threshold && relax_tye == 1)
                item.s = s * new_s * 2;
            end
            node  = node.next
        end
    end
end

function freq_inbreed(sel::Float64,freq::Float64,h::Float64,F::Float64) 
    #The fitnesses of the genotypes are A1A1 = 1; A1A2 = 1-h*(2*s); A2A2 = 1-(2*s), where A2 is the derived allele. Under any conditions that the program is run, the value of s for any segregating site cannot exceed 0.5.
    # parameters struct assume positive values of s to concieve positive selection. Changing to negative value
    s::Float64=-sel
    return (((1.0 - s) * (freq * freq + F * freq * (1.0 - freq))) + ((1.0 - h * s) * freq * (1.0 - freq) * (1.0 - F))) / (((1.0 - freq) * (1.0 - freq) + (1.0 - freq) * F * freq) + ((1.0 - h * s) * 2.0 * freq * (1.0 - freq) * (1.0 - F)) + ((1.0 - s) * (freq * freq + F * freq * (1.0 - freq))))
end

function add_mutation!(mutation_list::LinkedList{mutation},r::Ptr{gsl_rng},N::Float64,h::Float64,s::Float64,θ::Float64,freq::Float64,dfe::String,param_one::Float64,param_two::Float64,s_mult::Float64,n_anc::Int64,age::Int64,trajectories::Array{Int64,1},relax::Bool)

    num_mut::Int64 = Int(ran_poisson(r,θ/2.0))
    count_mut::Int64 = length(mutation_list)

    local s_value::Float64 = 0.0
    @inbounds for x::Int64=1:num_mut
        if dfe == "point"
            s_value = s;
        elseif dfe == "gamma"
            ##Use this for gamma distribution in Boyko 2008:
            gamma = ran_gamma(r, param_one, param_two * s_mult);
            # note, this is scale for a Nanc=1000, using the boyko params
            s_value = - gamma / (n_anc * 2);
        elseif dfe == "beta"
            gamma = ran_beta(r, param_one, param_two * s_mult);
            s_value = - gamma / (n_anc * 2);
        end

        count_mut += 1;
        mut        = mutation(frequency = freq, h = h, s = s_value*2, count_samp = 0.0, age = age, num = count_mut, type=1)
        push!(mutation_list,mut);
        
        if(!isempty(trajectories) && count_mut in trajectories)
            # @show count_mut
            # i::Int64 = findfirst(isequal(count_mut),trajectories);
            # @printf(io, "%d\t%lf\n", trajectories[i], (1.0/N));
            push!(trajectories_output[count_mut],(1.0/N))
        end
    end
end

function drift_sel!(mutation_list::LinkedList{mutation},r::Ptr{gsl_rng},N::Float64,F::Float64,trajectories::Array{Int64,1})

    l = 0;
    f = 0;

    if(!isempty(mutation_list))

        # Index to delete node
        node=mutation_list.node.next;
        while (node.next != mutation_list.node.next)
            item        = node.data;

            freq            = item.frequency;
            s::Float64      = item.s;
            h::Float64      = item.h;
            freq::Float64   = freq_inbreed(s,freq,h,F);
            
            num       = item.num;

            #count::Int64    = rand(Binomial(N,freq));
            count::Int64    = ran_binomial(r,freq,N)
            item.count_samp = count;
            freq            = count / N;
            item.frequency  = freq;

            if (freq > 0.0 && freq < 1.0)   
                if(!isempty(trajectories) && num in trajectories)
                    # @show num
                    # i = findfirst(isequal(num),trajectories);
                    # @printf(io, "%d\t%lf\n", trajectories[i], freq);
                    push!(trajectories_output[num],freq)
                end
                # Move node
                node = node.next
            else
                if(freq == 0)
                    l +=1;
                elseif(freq == 1)
                    f +=1;
                end
                # Delete node
                deleteat!(mutation_list,node);
                # Move node
                node = node.next
            end
        end
        return(l,f)
    else
        return(l,f)
    end
end

function burnin(param::recipe, r::Ptr{gsl_rng})

    @unpack N, θ, s, h, dfe, param_one, param_two, s_mult, n_anc, trajectories, relax = param

    n_size::Int64 = N[1] - 1
    s_size::Int64 = N[1] - 1
    sel::Float64 = s[1]

    theoretical_sfs = zeros(Float64, n_size, 2)
    number_of_mutations::Int64 = 0

    mutation_list_burnin = LinkedList{mutation}()

    @inbounds for j::Int64 = 2:n_size
        res, err = quadgk(x -> f_of_q_lambda(x, j, n_size, s_size, sel, θ * 2), 0.0, 1.0)
        age::Int64 = n_size * 10
        add_mutation!(mutation_list_burnin, r, Float64(n_size), h, sel, res, j / n_size, dfe, param_one, param_two, s_mult[1], n_anc, age, trajectories, relax)
    end

    return mutation_list_burnin
end

function f_of_q_lambda(x::Float64,j::Int64,N_burnin::Int64,sample_size::Int64,point_sel::Float64,theta::Float64)

    gamma = N_burnin * point_sel;

    if (abs(gamma) > 1.0e-7)
        
        sfs_function_term_one = (1-exp(-2*N_burnin*(-point_sel)*(1-x)))/(1-exp(-2*N_burnin*(-point_sel)));
        sfs_function_term_two::Float64 = (2/(x*(1-x)));

        binomial_value::Float64 = ran_binomial_pdf(j,x,sample_size);
        # binomial_value::Float64 = pdf(Binomial(sample_size,x), j);

        f = theta/2 * sfs_function_term_one * sfs_function_term_two * binomial_value;

    else
        f = theta/2 * 2/x * ran_binomial_pdf(j,x,sample_size);
    end
    return  f;
end
