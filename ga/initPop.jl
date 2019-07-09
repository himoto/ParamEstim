function getInitialPopulation(
    n_population::Int64,
    n_gene::Int64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64}
    )::Matrix{Float64}
    population::Matrix{Float64} = Inf.*ones(n_population,n_gene+1);
    println("initpop\n");

    for i = 1:n_population
        while isinf(population[i,end]) || isnan(population[i,end])
            for j = 1:n_gene
                population[i,j] = rand();
            end

            population[i,end] = getFitness(population[i,1:n_gene],searchIdx,searchRegion);

        end
        
        print(@sprintf("%d / %d\r", i, n_population));
        flush(stdout);
    end
    print(@sprintf("%d / %d\n",n_population,n_population));

    population = sortslices(population,dims=1,by=x->x[end]);

    return population
end