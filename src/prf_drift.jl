
function relax_selection!(
    mutation_list::Vector{Mutation},
    new_s::Float64,
    sel_threshold::Float64,
    relax_type::Int
)
    isempty(mutation_list) && return nothing

    @inbounds for i in eachindex(mutation_list)
        m = mutation_list[i]
        s = m.s

        if (s / 2) <= sel_threshold
            s_new = (relax_type == 0) ? (new_s * 2) :
                    (relax_type == 1) ? (s * new_s * 2) :
                    s

            if s_new != s
                mutation_list[i] = Mutation(m.frequency, m.h, s_new, m.count_samp, m.age, m.num, m.type)
            end
        end
    end

    return nothing
end



@inline function freq_inbreed(sel::Float64, p::Float64, h::Float32, F::Float64)
    # pq = p(1-p)
    pq  = p * (1.0 - p)

    # P22 homozygote for allele-2, P12 = half heterozygotes
    P22 = p*p + F * pq
    P12 = pq * (1.0 - F)

    # A = contribution of selected genotypes to allele-2 numerator term
    A   = P22 + h * P12

    # Numerator: p + sel*A  (because P22 + P12 == p)
    num = p + sel * A

    # Denominator: 1 + sel*(P22 + 2h*P12)  (using P11 = 1 - P22 - 2P12)
    den = 1.0 + sel * (P22 + 2.0*h*P12)

    # Optional safety clamp for numerical noise (keeps ran_binomial inputs valid)
    x = num / den
    return ifelse(x < 0.0, 0.0, ifelse(x > 1.0, 1.0, x))
end


@inline sample_s(dfe_dist::Float64,n_anc::Int64) = dfe_dist
@inline sample_s(dfe_dist::Vector{Float64},n_anc::Int64) = dfe_dist[1]
@inline sample_s(dfe_dist::Gamma{Float64}, n_anc::Int64) = -rand(dfe_dist) / (n_anc * 2)
@inline sample_s(dfe_dist::LogNormal{Float64}, n_anc::Int64) = -rand(dfe_dist) / (n_anc * 2)
@inline sample_s(dfe_dist::Beta{Float64}, n_anc::Int64) = -rand(dfe_dist)


function add_mutation!(
    mutation_list::Vector{Mutation},
    N_F::Float64,
    h::Float64,
    θ::Poisson,
    freq::Float64,
    dfe_dist,
    n_anc::Int,
    age::Int,
    state::Trajectories,
    r::Ptr{gsl_rng},
)
    num_mut = rand(θ) 
    count_mut = length(mutation_list)
    N_inv = 1.0 / N_F

    # Reduce reallocations during bursts of new mutations
    #sizehint!(mutation_list, count_mut + num_mut)

    @inbounds for _ = 1:num_mut
        s_value = sample_s(dfe_dist, n_anc)
        # C code cap s to 0.5
        # s_value = min(s_value, 0.5)
        
        state.next_mut_id += Int32(1)
        mut_id = state.next_mut_id

        # store 2s, consistent with your current struct design
        push!(mutation_list, Mutation(freq, Float32(h), s_value * 2, 0, UInt32(age), mut_id, 1))

        # trajectory tracking: if this ID is requested, append the starting freq
        v = get(state.trajectories_output, mut_id, nothing)
        if v !== nothing
            push!(v, freq)
        end

    end

    return nothing
end


function add_mutation!(
    mutation_list::Vector{Mutation},
    N_F::Float64,
    h::Float64,
    θ::Poisson,
    freq::Float64,
    dfe_dist,
    n_anc::Int,
    age::Int,
    state::Nothing,
    r::Ptr{gsl_rng},
)
    num_mut = rand(θ)
    count_mut = length(mutation_list)

    # Reduce reallocations during bursts of new mutations
    #sizehint!(mutation_list, count_mut + num_mut)

    @inbounds for _ = 1:num_mut
        s_value = sample_s(dfe_dist, n_anc)
        # C code cap s to 0.5
        # s_value = min(s_value, 0.5)
        
        count_mut += 1

        push!(mutation_list, Mutation(freq, h, s_value * 2, 0, age, count_mut, 1))

    end

    return nothing
end

function drift_sel!(
    mutation_list::Vector{Mutation},
    r::Ptr{gsl_rng},
    N::Float64,
    F::Float64,
    state::Trajectories
)
    l = 0
    f = 0
    isempty(mutation_list) && return (l, f)

    N_inv = 1.0 / N
    i = 1

    @inbounds while i <= length(mutation_list)
        m = mutation_list[i]

        freq = freq_inbreed(m.s, m.frequency, m.h, F)

        count = ran_binomial(r, freq, N)
        freq  = count * N_inv

        if 0.0 < freq < 1.0
            # update in place by replacement (isbits copy)
            mutation_list[i] = Mutation(freq, m.h, m.s, count, m.age, m.num, m.type)


            v = get(state.trajectories_output, m.num, nothing)
            if v !== nothing
                push!(v, freq)
            end

            i += 1
        else
            l += (freq == 0.0)
            f += (freq == 1.0)

            # swap-delete: overwrite current with last, then pop
            mutation_list[i] = mutation_list[end]
            pop!(mutation_list)

            # do not increment i; process the swapped-in element next
        end
    end

    return (l, f)
end

function drift_sel!(
    mutation_list::Vector{Mutation},
    r::Ptr{gsl_rng},
    N::Float64,
    F::Float64,
    state::Nothing
)
    l = 0
    f = 0
    isempty(mutation_list) && return (l, f)

    N_inv = 1.0 / N
    i = 1

    @inbounds while i <= length(mutation_list)
        m = mutation_list[i]

        freq = freq_inbreed(m.s, m.frequency, m.h, F)

        count = ran_binomial(r, freq, N)
        freq  = count * N_inv

        if 0.0 < freq < 1.0
            # update in place by replacement (isbits copy)
            mutation_list[i] = Mutation(freq, m.h, m.s, count, m.age, m.num, m.type)
            i += 1
        else
            l += (freq == 0.0)
            f += (freq == 1.0)

            # swap-delete: overwrite current with last, then pop
            mutation_list[i] = mutation_list[end]
            pop!(mutation_list)

            # do not increment i; process the swapped-in element next
        end
    end

    return (l, f)
end


function burnin(param::recipe,dfe_dist::T) where T <: Union{Vector{Float64}, Gamma{Float64}, Beta{Float64}}

    @unpack N, θ, s, h, dfe, param_one, param_two, s_mult, n_anc, trajectories, relax = param

    n_size::Int64 = N[1] - 1
    s_size::Int64 = N[1] - 1
    sel::Float64 = s[1]

    theoretical_sfs = zeros(Float64, n_size, 2)
    number_of_mutations::Int64 = 0

    mutation_list_burnin = LinkedList{mutation}()

    @inbounds for j::Int64 = 2:n_size
        res, err = quadgk(x -> f_of_q_lambda(x, j, n_size, s_size, sel, θ), 0.0, 1.0)
        age::Int64 = -1
        # add_mutation!(mutation_list_burnin, Float64(n_size), h, Poisson(res), j / n_size, dfe_dist, n_anc, age, relax)
        add_mutation!(mutation_list,N_F,h,θ_dist,N_freq,dfe_dist,n_anc,age,nothing,r)
    end

    return mutation_list_burnin
end

function f_of_q_lambda(x::Float64,j::Int64,N_burnin::Int64,sample_size::Int64,point_sel::Float64,theta::Float64)

    gamma = N_burnin * point_sel;

    if (abs(gamma) > 1.0e-7)
        
        sfs_function_term_one = (1-exp(-2*N_burnin*(-point_sel)*(1-x)))/(1-exp(-2*N_burnin*(-point_sel)));
        sfs_function_term_two::Float64 = (2/(x*(1-x)));

        # binomial_value::Float64 = ran_binomial_pdf(j,x,sample_size);
        binomial_value::Float64 = pdf(Binomial(sample_size,x), j);

        f = theta/2 * sfs_function_term_one * sfs_function_term_two * binomial_value;

    else
        # f = theta/2 * 2/x * ran_binomial_pdf(j,x,sample_size);
        f = theta/2 * 2/x * pdf(Binomial(sample_size,x), j);
    end
    return  f;
end
