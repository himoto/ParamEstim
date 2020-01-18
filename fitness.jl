# Residual Sum of Squares
function compute_objval_rss(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64
    error::Float64 = 0.0

    for i in eachindex(exp_data)
        error += (sim_data[i]-exp_data[i])^2
    end

    return error
end


# Cosine similarity
function compute_objval_cos(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64

    error::Float64 = 1.0 - dot(sim_data,exp_data)/(norm(sim_data)*norm(exp_data))

    return error
end


# Define an objective function to be minimized.
function objective(individual_gene::Vector{Float64}, search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                    search_region::Matrix{Float64})::Float64
    p,u0 = update_param(
        individual_gene::Vector{Float64},
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
        search_region::Matrix{Float64}
    )

    # constraints --------------------------------------------------------------
    p[C.V6] = p[C.V5]
    p[C.Km6] = p[C.Km5]
    p[C.KimpDUSP] = p[C.KimDUSP]
    p[C.KexpDUSP] = p[C.KexDUSP]
    p[C.KimpcFOS] = p[C.KimFOS]
    p[C.KexpcFOS] = p[C.KexFOS]
    p[C.p52] = p[C.p47]
    p[C.m52] = p[C.m47]
    p[C.p53] = p[C.p48]
    p[C.p54] = p[C.p49]
    p[C.m54] = p[C.m49]
    p[C.p55] = p[C.p50]
    p[C.p56] = p[C.p51]
    p[C.m56] = p[C.m51]
    # --------------------------------------------------------------------------

    if Sim.simulate!(p,u0) isa Nothing
        error::Vector{Float64} = zeros(length(observables))
        for (i,target) in enumerate(observables)
            exp_t::Vector{Float64} = Exp.get_timepoint(i)
            norm_max::Float64 = maximum(Sim.simulations[obs_idx(target),:,:])
            if isassigned(Exp.experiments,obs_idx(target))
                error[i] = compute_objval_rss(
                    diff_sim_and_exp(
                        Sim.simulations[obs_idx(target),:,:],
                        Exp.experiments[obs_idx(target)],
                        exp_t,
                        Sim.conditions,
                        sim_norm_max=norm_max,
                        exp_norm_max=1.0
                    )...
                )
            end
        end
    
        return sum(error)
    else
        return Inf
    end
end