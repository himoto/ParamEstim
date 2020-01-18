const observables = [
    "Phosphorylated_MEKc"
    "Phosphorylated_ERKc"
    "Phosphorylated_RSKw"
    "Phosphorylated_CREBw"
    "dusp_mRNA"
    "cfos_mRNA"
    "cFos_Protein"
    "Phosphorylated_cFos"
];

function obs_idx(observable_name::String)::Int

    return findfirst(isequal(observable_name),observables)
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