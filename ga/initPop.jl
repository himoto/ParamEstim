function getInitialPopulation(n_population::Int64, n_gene::Int64, searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
                                searchRegion::Matrix{Float64})::Matrix{Float64}
    population::Matrix{Float64} = Inf.*ones(n_population,n_gene+1)
    
    println("Generating the initial population...\n")
    flush(stdout)

    for i = 1:n_population
        while isinf(population[i,end]) || isnan(population[i,end])
            for j = 1:n_gene
                population[i,j] = rand()
            end
            population[i,end] = objective(
                population[i,1:n_gene],searchIdx,searchRegion
            )
        end
        print(@sprintf("%d / %d\n", i, n_population))
        flush(stdout)
    end
    print("\n######\n\n")

    population = sortslices(population,dims=1,by=x->x[end])

    return population
end


function getInitialPopulation_continue(nthParamSet::Int64, n_population::Int64, n_gene::Int64,
                                        searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
                                        searchRegion::Matrix{Float64}, p0_bounds::Vector{Float64})::Matrix{Float64}
    generation::Int64 = readdlm(
        "./fitparam/$nthParamSet/generation.dat"
    )[1,1]
    bestIndiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat",nthParamSet,generation
        )
    )[:,1]
    population::Matrix{Float64} = Inf.*ones(n_population,n_gene+1)

    println("Generating the initial population...\n")
    flush(stdout)

    for i = 1:n_population
        while isinf(population[i,end]) || isnan(population[i,end])
            for j = 1:n_gene
                population[i,j] = (
                    log10(
                        bestIndiv[j]*10^(
                            rand()*log10(p0_bounds[2]/p0_bounds[1])+log10(p0_bounds[1])
                        )
                    ) - searchRegion[1,j]
                )/(searchRegion[2,j] - searchRegion[1,j])
            end
            for j=1:n_gene
                if population[i,j] > 1.0
                    population[i,j] = 1.0
                elseif population[i,j] < 0.0
                    population[i,j] = 0.0
                end
            end
            population[i,end] = objective(population[i,1:n_gene],searchIdx,searchRegion)
        end
        print(@sprintf("%d / %d\n", i, n_population))
        flush(stdout)
    end
    print("\n----------------------------------------\n\n")

    population = sortslices(population,dims=1,by=x->x[end])

    return population
end