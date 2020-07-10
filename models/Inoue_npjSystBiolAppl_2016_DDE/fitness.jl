# Residual Sum of Squares
function compute_objval_rss(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64
    error::Float64 = 0.0

    @simd for i in eachindex(exp_data)
        @inbounds error += (sim_data[i] - exp_data[i])^2
    end

    return error
end


# Cosine similarity
function compute_objval_cos(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64

    error::Float64 = 1.0 - dot(sim_data,exp_data)/(norm(sim_data)*norm(exp_data))

    return error
end


function conditions_index(condition_name::String)::Int

    return findfirst(isequal(condition_name),Sim.conditions)
end


function diff_sim_and_exp(
        sim_matrix::Matrix{Float64},
        exp_dict::Dict{String,Array{Float64,1}},
        exp_timepoint::Vector{Float64},
        conditions::Vector{String};
        sim_norm_max::Float64)
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (idx,condition) in enumerate(conditions)
        if condition in keys(exp_dict)
            append!(sim_result,sim_matrix[Int.(exp_timepoint.+1),idx])
            append!(exp_result,exp_dict[condition])
        end
    end

    return (sim_result./sim_norm_max, exp_result)
end


# Define an objective function to be minimized.
function objective(indiv_gene::Vector{Float64})::Float64
    indiv::Vector{Float64} = decode_gene2val(indiv_gene)

    (p,u0) = update_param(indiv)

    if Sim.simulate!(p,u0) isa Nothing
        error::Vector{Float64} = zeros(length(observables))
        for (i,target) in enumerate(observables)
            if isassigned(Exp.experiments,observables_index(target))
                error[i] = compute_objval_cos(
                    diff_sim_and_exp(
                        Sim.simulations[observables_index(target),:,:],
                        Exp.experiments[observables_index(target)],
                        Exp.get_timepoint(i),
                        Sim.conditions,
                        sim_norm_max=ifelse(
                            Sim.normalization,
                            maximum(Sim.simulations[observables_index(target),:,:]),
                            1.0
                        ),
                    )...
                )
            end
        end
        return sum(error)
    else
        return Inf
    end
end