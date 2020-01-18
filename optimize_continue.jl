include("ParamEstim.jl")
using .ParamEstim

function optimize_continue(nth_param_set::Int64)

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = search_parameter_index()
    search_region::Matrix{Float64} = get_search_region()

    max_generation::Int64 = 10000
    n_population::Int64 = 5*size(search_region,2)
    n_children::Int64 = 50
    n_gene::Int64 = size(search_region,2)
    allowable_error::Float64 = 0.35

    p0_bounds = [0.1,10.0]  # [lower_bound,upper_bound]

    if !isdir("./fitparam/$nth_param_set")
        mkdir("./fitparam/$nth_param_set")
        
        (bestIndiv,bestFitness) = ga_v2(
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            search_idx,
            search_region
        )
    else
        (bestIndiv,bestFitness) = ga_v2_continue(
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            search_idx,
            search_region,
            p0_bounds
        )
    end
end

###
optimize_continue(parse(Int64,ARGS[1]))