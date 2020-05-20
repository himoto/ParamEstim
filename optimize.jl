include("ParamEstim.jl")
using .ParamEstim

function optimize(nth_param_set::Int64)
    if !isdir("./fitparam")
        mkdir("./fitparam")
    end

    try
        files = readdir(
            "./fitparam/$nth_param_set"
        )
        for file in files
            if occursin(".dat",file)
                rm(
                    "./fitparam/$nth_param_set/$file"
                )
            end
        end
    catch
        mkdir("./fitparam/$nth_param_set")
    end

    search_region::Matrix{Float64} = get_search_region()

    max_generation::Int64 = 10000
    n_population::Int64 = 15*size(search_region, 2)
    n_children::Int64 = 50
    n_gene::Int64 = size(search_region, 2)
    allowable_error::Float64 = 0.0

    (best_indiv, best_fitness) = ga_v2(
        nth_param_set,
        max_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error,
        search_region
    )
end


if abspath(PROGRAM_FILE) == @__FILE__
    optimize(parse(Int64,ARGS[1]))
end