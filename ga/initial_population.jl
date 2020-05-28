function get_initial_population(n_population::Int64, n_gene::Int64,
                                search_rgn::Matrix{Float64})::Matrix{Float64}
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
                population[i,1:n_gene], search_rgn
            )
        end
        print(
            @sprintf(
                "%d / %d\n", i, n_population
            )
        )
        flush(stdout)
    end
    print(
        "\n----------------------------------------\n\n"
    )
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end


function encode_bestindiv2randgene(j::Int64, best_indiv::Vector{Float64}, 
                                    search_rgn::Matrix{Float64},
                                    p0_bounds::Vector{Float64})::Float64
    rand_gene::Float64 = (
        log10(
            best_indiv[j]*10^(
                rand()*log10(p0_bounds[2]/p0_bounds[1]) + log10(p0_bounds[1])
            )
        ) - search_rgn[1,j]
    ) / (
        search_rgn[2,j] - search_rgn[1,j]
    )

    return rand_gene
end


function get_initial_population_continue(nthParamSet::Int64, n_population::Int64, n_gene::Int64,
                                            search_rgn::Matrix{Float64},
                                            p0_bounds::Vector{Float64})::Matrix{Float64}
    generation::Int64 = readdlm(
        "./fitparam/$nthParamSet/generation.dat"
    )[1,1]
    best_indiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat", nthParamSet, generation
        )
    )[:,1]
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
                population[i,j] = encode_bestindiv2randgene(
                    j, best_indiv, search_rgn, p0_bounds
                )
                if population[i,j] > 1.0
                    population[i,j] = 1.0
                elseif population[i,j] < 0.0
                    population[i,j] = 0.0
                end
            end
            population[i,end] = objective(
                population[i,1:n_gene], search_rgn
            )
        end
        print(
            @sprintf(
                "%d / %d\n", i, n_population
            )
        )
        flush(stdout)
    end
    print(
        "\n----------------------------------------\n\n"
    )
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end