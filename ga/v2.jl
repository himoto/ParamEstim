function ga_v2(
        nth_param_set::Int64,
        max_generation::Int64,
        n_population::Int64,
        n_children::Int64,
        n_gene::Int64,
        allowable_error::Float64)::Tuple{Array{Float64, 1}, Float64}
    if n_population < n_gene + 2
        error("n_population must be larger than $(n_gene+2)")
    end

    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(3*n_population)

    population = get_initial_population(n_population, n_gene)
    N0[1] = population[1, end]
    println(
        @sprintf(
            "Generation%d: Best Fitness = %.6e", 1, population[1, end]
        )
    )
    flush(stdout)

    best_indiv = decode_gene2val(population[1, 1:n_gene])
    best_fitness = population[1, end]

    f = open("./fitparam/$nth_param_set/fit_param1.dat", "w")
    for val in best_indiv
        write(f, @sprintf("%.6e\n", val))
    end
    close(f)

    open("./fitparam/$nth_param_set/generation.dat", "w") do f
        write(f, @sprintf("%d", 1))
    end

    open("./fitparam/$nth_param_set/best_fitness.dat", "w") do f
        write(f, @sprintf("%.6e", best_fitness))
    end

    if population[1, end] <= allowable_error
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]

        return best_indiv, best_fitness
    end

    generation::Int64 = 2
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        population = converging!(ip, population, n_population, n_gene)
        population = local_search!(ip, population, n_population, n_children, n_gene)
        if N_iter > 1
            for _ in 1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                population = converging!(ip, population, n_population, n_gene)
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1, end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        else
            N0[generation%length(N0)] = population[1, end]
        end

        println(
            @sprintf(
                "Generation%d: Best Fitness = %.6e", generation, population[1, end]
            )
        )
        flush(stdout)
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        if population[1, end] < best_fitness
            f = open(
                "./fitparam/$nth_param_set/fit_param$generation.dat", "w"
            )
            for val in best_indiv
                write(f, @sprintf("%.6e\n", val))
            end
            close(f)

            open("./fitparam/$nth_param_set/generation.dat", "w") do f
                write(f, @sprintf("%d", generation))
            end
        end
        best_fitness = population[1, end]

        open("./fitparam/$nth_param_set/best_fitness.dat", "w") do f
            write(f, @sprintf("%.6e", best_fitness))
        end

        if population[1, end] <= allowable_error
            best_indiv = decode_gene2val(population[1, 1:n_gene])
            best_fitness = population[1, end]

            return best_indiv, best_fitness
        end

        open("./fitparam/$nth_param_set/count_num.dat", "w") do f
            write(f, @sprintf("%d", generation))
        end

        generation += 1
    end

    best_indiv = decode_gene2val(population[1, 1:n_gene])
    best_fitness = population[1, end]

    return best_indiv, best_fitness
end


function ga_v2_continue(
        nth_param_set::Int64,
        max_generation::Int64,
        n_population::Int64,
        n_children::Int64,
        n_gene::Int64,
        allowable_error::Float64, 
        p0_bounds::Vector{Float64})::Tuple{Array{Float64, 1}, Float64}
    if n_population < n_gene + 2
        error("n_population must be larger than $(n_gene+2)")
    end
    
    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(3*n_population)

    count::Int64 = readdlm(
        "./fitparam/$nth_param_set/count_num.dat"
    )[1, 1]
    best_generation::Int64 = readdlm(
        "./fitparam/$nth_param_set/generation.dat"
    )[1, 1]
    best_indiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat", nth_param_set, best_generation
        )
    )[:, 1]
    best_indiv_gene::Vector{Float64} = encode_val2gene(best_indiv)
    best_fitness::Float64 = objective(best_indiv_gene)

    population = get_initial_population_continue(
        nth_param_set, n_population, n_gene, p0_bounds
    )
    if best_fitness < population[1, end]
        for i in 1:n_gene
            @inbounds population[1, i] = best_indiv_gene[i]
        end
        population[1, end] = best_fitness
    else
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]
        open("./fitparam/$nth_param_set/fit_param$count.dat", "w") do f
            for i=1:n_gene
                write(f, @sprintf("%.6e", best_indiv[i]))
            end
        end
    end

    N0[1] = population[1, end]

    println(
        @sprintf(
            "Generation%d: Best Fitness = %.6e", count + 1, population[1, end]
        )
    )
    flush(stdout)

    if population[1, end] <= allowable_error
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]

        return best_indiv, best_fitness
    end

    generation::Int64 = 2 + count
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        population = converging!(ip, population, n_population, n_gene)
        population = local_search!(ip, population, n_population, n_children, n_gene)
        if N_iter > 1
            for _ in 1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                population = converging!(ip, population, n_population, n_gene)
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1, end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        else
            N0[generation%length(N0)] = population[1, end]
        end

        println(
            @sprintf(
                "Generation%d: Best Fitness = %.6e", generation, population[1, end]
            )
        )
        flush(stdout)
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        if population[1, end] < best_fitness
            f = open(
                @sprintf(
                    "./fitparam/%d/fit_param%d.dat", nth_param_set, generation
                ), "w"
            )
            for val in best_indiv
                write(f, @sprintf("%.6e\n", val))
            end
            close(f)

            open("./fitparam/$nth_param_set/generation.dat", "w") do f
                write(f, @sprintf("%d", generation))
            end
        end
        best_fitness = population[1, end]

        open("./fitparam/$nth_param_set/best_fitness.dat", "w") do f
            write(f, @sprintf("%.6e", best_fitness))
        end

        if population[1, end] <= allowable_error
            best_indiv = decode_gene2val(population[1, 1:n_gene])
            best_fitness = population[1, end]

            return best_indiv, best_fitness
        end

        open("./fitparam/$nth_param_set/count_num.dat", "w") do f
            write(f, @sprintf("%d", generation))
        end

        generation += 1
    end

    best_indiv = decode_gene2val(population[1, 1:n_gene])
    best_fitness = population[1, end]

    return best_indiv, best_fitness
end