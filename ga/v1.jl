function gaV1(
    nthParamSet::Int64,
    n_generation::Int64,
    n_population::Int64,
    n_children::Int64,
    n_gene::Int64,
    allowable_error::Float64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64}
    )::Tuple{Array{Float64,1},Float64}
    population = getInitialPopulation(n_population,n_gene,searchIdx,searchRegion);
    print(@sprintf("Generation%d: Best Fitness = %e\n",1,population[1,end]));
    flush(stdout);

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    f = open("../FitParam/$nthParamSet/fitParam1.dat", "w");
    for i=1:length(searchIdx[1])
        write(f,@sprintf("%e\n",bestIndiv[i]));
    end
    for i=1:length(searchIdx[2])
        write(f,@sprintf("%e\n",bestIndiv[i+length(searchIdx[1])]));
    end
    close(f);

    open("../FitParam/$nthParamSet/generation.dat", "w") do f
        write(f,@sprintf("%d",1));
    end

    open("../FitParam/$nthParamSet/bestFitness.dat", "w") do f
        write(f,@sprintf("%e",bestFitness));
    end

    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        bestFitness = population[1,end];
        return bestIndiv,bestFitness
    end

    for i = 2:n_generation
        population = mggVariant(population,n_population,n_children,n_gene,searchIdx,searchRegion);
        print(@sprintf("Generation%d: Best Fitness = %e\n",i,population[1,end]));
        flush(stdout);
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);

        if population[1,end] < bestFitness
            f = open("../FitParam/$nthParamSet/fitParam$i.dat", "w");
            for i=1:length(searchIdx[1])
                write(f,@sprintf("%e\n",bestIndiv[i]));
            end
            for i=1:length(searchIdx[2])
                write(f,@sprintf("%e\n",bestIndiv[i+length(searchIdx[1])]));
            end
            close(f);

            open("../FitParam/$nthParamSet/generation.dat", "w") do f
                write(f,@sprintf("%d",i));
            end
        end

        bestFitness = population[1,end];

        open("../FitParam/$nthParamSet/bestFitness.dat", "w") do f
            write(f,@sprintf("%e",bestFitness));
        end

        if population[1,end] <= allowable_error
            bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
            bestFitness = population[1,end];
            return bestIndiv,bestFitness
        end
    end

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    return bestIndiv,bestFitness
end