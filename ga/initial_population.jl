function get_initial_population(
        n_population::Int64,
        n_gene::Int64)::Matrix{Float64}
    println(
        "Generating the initial population...\n"
    )
    flush(stdout)
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    for i = 1:n_population
        while !isfinite(population[i,end])
            for j = 1:n_gene
                population[i,j] = rand()
            end
            population[i,end] = objective(
                population[i,1:n_gene]
            )
        end
        println(
            @sprintf(
                "%d / %d", i, n_population
            )
        )
        flush(stdout)
    end
    println(
        "\n----------------------------------------\n"
    )
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end


function get_initial_population_continue(
        nth_param_set::Int64,
        n_population::Int64,
        n_gene::Int64,
        p0_bounds::Vector{Float64})::Matrix{Float64}
    generation::Int64 = readdlm(
        "./fitparam/$nth_param_set/generation.dat"
    )[1,1]
    best_indiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat", nth_param_set, generation
        )
    )[:,1]
    println(
        "\n----------------------------------------\n"*
        "Generating the initial population...\n"
    )
    flush(stdout)
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    for i = 1:n_population
        while !isfinite(population[i,end])
            for j = 1:n_gene
                population[i,j] = encode_bestIndivVal2randGene(
                    j, best_indiv, p0_bounds
                )
                if population[i,j] > 1.0
                    population[i,j] = 1.0
                elseif population[i,j] < 0.0
                    population[i,j] = 0.0
                end
            end
            population[i,end] = objective(
                population[i,1:n_gene]
            )
        end
        println(
            @sprintf(
                "%d / %d", i, n_population
            )
        )
        flush(stdout)
    end
    println(
        "\n----------------------------------------\n"
    )
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end