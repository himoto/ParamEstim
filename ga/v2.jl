function gaV2(nthParamSet::Int64, max_generation::Int64, n_population::Int64, n_children::Int64,
                n_gene::Int64, allowable_error::Float64, searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
                searchRegion::Matrix{Float64})::Tuple{Array{Float64,1},Float64}
    if n_population < n_gene+2
        error("n_population must be larger than $(n_gene+2)")
    end

    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(2*n_population)

    population = getInitialPopulation(
        n_population,n_gene,searchIdx,searchRegion
    )
    N0[1] = population[1,end]
    print(
        @sprintf(
            "Generation%d: Best Fitness = %.6e\n",1,population[1,end]
        )
    )
    flush(stdout)

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
    bestFitness = population[1,end]

    f = open("./fitparam/$nthParamSet/fit_param1.dat", "w")
    for i in eachindex(searchIdx[1])
        write(f,@sprintf("%.6e\n",bestIndiv[i]))
    end
    for i in eachindex(searchIdx[2])
        write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]))
    end
    close(f)

    open("./fitparam/$nthParamSet/generation.dat", "w") do f
        write(f,@sprintf("%d",1))
    end

    open("./fitparam/$nthParamSet/best_fitness.dat", "w") do f
        write(f,@sprintf("%.6e",bestFitness))
    end

    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
        bestFitness = population[1,end]
        return bestIndiv,bestFitness
    end

    generation::Int64 = 2
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        ip, population = converging(
            ip,population,n_population,n_gene,searchIdx,searchRegion
        )
        ip, population = localsearch(
            ip,population,n_population,n_children,n_gene,searchIdx,searchRegion
        )
        if N_iter > 1
            for j=1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                ip, population = converging(
                    ip,population,n_population,n_gene,searchIdx,searchRegion
                )
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1,end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        elseif generation%length(N0) == 1
            N0 = zeros(2*n_population)
            N0[1] = population[1,end]
        else
            N0[generation%length(N0)] = population[1,end]
        end

        print(
            @sprintf(
                "Generation%d: Best Fitness = %.6e\n",generation,population[1,end]
            )
        )
        flush(stdout)
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
        if population[1,end] < bestFitness
            f = open("./fitparam/$nthParamSet/fit_param$generation.dat", "w")
            for i in eachindex(searchIdx[1])
                write(f,@sprintf("%.6e\n",bestIndiv[i]))
            end
            for i in eachindex(searchIdx[2])
                write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]))
            end
            close(f)

            open("./fitparam/$nthParamSet/generation.dat", "w") do f
                write(f,@sprintf("%d",generation))
            end
        end
        bestFitness = population[1,end]

        open("./fitparam/$nthParamSet/best_fitness.dat", "w") do f
            write(f,@sprintf("%.6e",bestFitness))
        end

        if population[1,end] <= allowable_error
            bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
            bestFitness = population[1,end]
            return bestIndiv,bestFitness
        end

        open("./fitparam/$nthParamSet/count_num.dat", "w") do f
            write(f,@sprintf("%d",generation))
        end
        generation += 1
    end
    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
    bestFitness = population[1,end]

    return bestIndiv,bestFitness
end


function gaV2_continue(nthParamSet::Int64, max_generation::Int64, n_population::Int64,
                        n_children::Int64, n_gene::Int64, allowable_error::Float64,
                        searchIdx::Tuple{Array{Int64,1},Array{Int64,1}}, searchRegion::Matrix{Float64},
                        p0_bounds::Vector{Float64})::Tuple{Array{Float64,1},Float64}
    if n_population < n_gene+2
        error("n_population must be larger than $(n_gene+2)")
    end
    
    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(2*n_population)

    count::Int64 = readdlm(
        "./fitparam/$nthParamSet/count_num.dat"
    )[1,1]
    bestGeneration::Int64 = readdlm(
        "./fitparam/$nthParamSet/generation.dat"
    )[1,1]
    bestIndiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat",nthParamSet,bestGeneration
        )
    )[:,1]
    bestFitness::Float64 = objective(
        (log10.(bestIndiv) .- searchRegion[1,:])./(searchRegion[2,:] .- searchRegion[1,:]),
        searchIdx,searchRegion
    )

    population = getInitialPopulation_continue(
        nthParamSet,n_population,n_gene,searchIdx,searchRegion,p0_bounds
    )
    if bestFitness < population[1,end]
        for i=1:n_gene
            population[1,i] = (
                log10(bestIndiv[i])-searchRegion[1,i])/(searchRegion[2,i]-searchRegion[1,i]
            )
        end
        population[1,end] = bestFitness
    else
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
        bestFitness = population[1,end]
        open("./fitparam/$nthParamSet/fit_param$count.dat", "w") do f
            for i=1:n_gene
                write(f,@sprintf("%.6e",bestIndiv[i]))
            end
        end
    end

    N0[1] = population[1,end]

    print(
        @sprintf(
            "Generation%d: Best Fitness = %.6e\n",count + 1,population[1,end]
        )
    )
    flush(stdout)

    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
        bestFitness = population[1,end]
        return bestIndiv,bestFitness
    end

    generation::Int64 = 2
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        ip, population = converging(
            ip,population,n_population,n_gene,searchIdx,searchRegion
        )
        ip, population = localsearch(
            ip,population,n_population,n_children,n_gene,searchIdx,searchRegion
        )
        if N_iter > 1
            for j=1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                ip, population = converging(
                    ip,population,n_population,n_gene,searchIdx,searchRegion
                )
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1,end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        elseif generation%length(N0) == 1
            N0 = zeros(2*n_population)
            N0[1] = population[1,end]
        else
            N0[generation%length(N0)] = population[1,end]
        end

        print(
            @sprintf(
                "Generation%d: Best Fitness = %.6e\n",generation + count,population[1,end]
            )
        )
        flush(stdout)
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)

        if population[1,end] < bestFitness
            f = open(
                @sprintf("./fitparam/%d/fit_param%d.dat",nthParamSet,generation + count), "w"
            )
            for i in eachindex(searchIdx[1])
                write(f,@sprintf("%.6e\n",bestIndiv[i]))
            end
            for i in eachindex(searchIdx[2])
                write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]))
            end
            close(f)

            open("./fitparam/$nthParamSet/generation.dat", "w") do f
                write(f,@sprintf("%d",generation + count))
            end
        end
        bestFitness = population[1,end]

        open("./fitparam/$nthParamSet/best_fitness.dat", "w") do f
            write(f,@sprintf("%.6e",bestFitness))
        end

        if population[1,end] <= allowable_error
            bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
            bestFitness = population[1,end]
            return bestIndiv,bestFitness
        end

        open("./fitparam/$nthParamSet/count_num.dat", "w") do f
            write(f,@sprintf("%d",generation + count))
        end
        generation += 1
    end
    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion)
    bestFitness = population[1,end]

    return bestIndiv,bestFitness
end