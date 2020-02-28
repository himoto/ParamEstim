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


function cond_idx(condition_name::String)::Int

    return findfirst(isequal(condition_name),Sim.conditions)
end


function diff_sim_and_exp(sim_matrix::Matrix{Float64},exp_dict::Dict{String,Array{Float64,1}},
                            exp_timepoint::Vector{Float64},conditions::Vector{String};
                            sim_norm_max::Float64,exp_norm_max::Float64)
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (i,condition) in enumerate(conditions)
        append!(sim_result,sim_matrix[Int.(exp_timepoint.+1),i])
        append!(exp_result,exp_dict[condition])
    end

    return (sim_result./sim_norm_max, exp_result./exp_norm_max)
end


# Define an objective function to be minimized.
function objective(individual_gene::Vector{Float64}, search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                    search_region::Matrix{Float64})::Float64
    p::Vector{Float64} = f_params()
    u0::Vector{Float64} = initial_values()
    
    indiv::Vector{Float64} = decode_gene2variable(individual_gene,search_region)

    for (i,j) in enumerate(search_idx[1])
        @inbounds p[j] = indiv[i]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds u0[j] = indiv[i+length(search_idx[1])]
    end

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
        #=
        error::Vector{Float64} = zeros(16)

        sim_norm_max = maximum(Sim.simulations[obs_idx("Phosphorylated_MEKc"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("Phosphorylated_MEKc"))
        error[1] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_MEKc"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_MEKc")]["EGF"]
        )
        error[2] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_MEKc"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_MEKc")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("Phosphorylated_ERKc"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("Phosphorylated_ERKc"))
        error[3] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_ERKc"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_ERKc")]["EGF"]
        )
        error[4] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_ERKc"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_ERKc")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("Phosphorylated_RSKw"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("Phosphorylated_RSKw"))
        error[5] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_RSKw"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_RSKw")]["EGF"]
        )
        error[6] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_RSKw"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_RSKw")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("Phosphorylated_CREBw"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("Phosphorylated_CREBw"))
        error[7] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_CREBw"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_CREBw")]["EGF"]
        )
        error[8] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_CREBw"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_CREBw")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("dusp_mRNA"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("dusp_mRNA"))
        error[9] = compute_objval_rss(
            Sim.simulations[obs_idx("dusp_mRNA"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("dusp_mRNA")]["EGF"]
        )
        error[10] = compute_objval_rss(
            Sim.simulations[obs_idx("dusp_mRNA"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("dusp_mRNA")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("cfos_mRNA"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("cfos_mRNA"))
        error[11] = compute_objval_rss(
            Sim.simulations[obs_idx("cfos_mRNA"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("cfos_mRNA")]["EGF"]
        )
        error[12] = compute_objval_rss(
            Sim.simulations[obs_idx("cfos_mRNA"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("cfos_mRNA")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("cFos_Protein"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("cFos_Protein"))
        error[13] = compute_objval_rss(
            Sim.simulations[obs_idx("cFos_Protein"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("cFos_Protein")]["EGF"]
        )
        error[14] = compute_objval_rss(
            Sim.simulations[obs_idx("cFos_Protein"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("cFos_Protein")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[obs_idx("Phosphorylated_cFos"),:,:])
        exp_t = Exp.get_timepoint(obs_idx("Phosphorylated_cFos"))
        error[15] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_cFos"), Int.(exp_t.+1), cond_idx("EGF")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_cFos")]["EGF"]
        )
        error[16] = compute_objval_rss(
            Sim.simulations[obs_idx("Phosphorylated_cFos"), Int.(exp_t.+1), cond_idx("HRG")]./sim_norm_max,
            Exp.experiments[obs_idx("Phosphorylated_cFos")]["HRG"]
        )
        =#
        return sum(error)
    else
        return Inf
    end
end