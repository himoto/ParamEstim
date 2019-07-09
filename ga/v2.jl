function gaV2(
    n_generation::Int64,
    n_population::Int64,
    n_children::Int64,
    n_gene::Int64,
    allowable_error::Float64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64}
    )::Tuple{Array{Float64,1},Float64}
    N_iter::Int64 = 1;
    N0::Vector{Float64} = zeros(2*n_population);

    population = getInitialPopulation(n_population,n_gene,searchIdx,searchRegion);
    N0[1] = population[1,end]
    print(@sprintf("Generation%d: Best Fitness = %e\n",1,population[1,end]));
    flush(stdout);

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    f = open("./FitParam/fitParam1.dat", "w");
    for i=1:length(searchIdx[1])
        write(f,@sprintf("%e\n",bestIndiv[i]));
    end
    for i=1:length(searchIdx[2])
        write(f,@sprintf("%e\n",bestIndiv[i+length(searchIdx[1])]));
    end
    close(f);

    open("./FitParam/generation.dat", "w") do f
        write(f,@sprintf("%d",1));
    end
    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        bestFitness = population[1,end];
        return bestIndiv,bestFitness
    end

    for i = 2:n_generation
        ip = randperm(n_population)[1:n_gene+2]
        ip, population = converging(ip,population,n_population,n_gene,searchIdx,searchRegion);
        ip, population = localsearch(ip,population,n_population,n_children,n_gene,searchIdx,searchRegion);
        if N_iter > 1
            for j=1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                ip, population = converging(ip,population,n_population,n_gene,searchIdx,searchRegion);
            end
        end

        if i%length(N0) == 0
            N0[end] = population[1,end];
            if N0[1] == N0[end]
                N_iter *= 2;
            else
                N_iter = 1;
            end
        elseif i%length(N0) == 1
            N0 = zeros(2*n_population);
            N0[1] = population[1,end];
        else
            N0[i%length(N0)] = population[1,end];
        end

        print(@sprintf("Generation%d: Best Fitness = %e\n",i,population[1,end]));
        flush(stdout);
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        if population[1,end] < bestFitness
            f = open("./FitParam/fitParam$i.dat", "w");
            for i=1:length(searchIdx[1])
                write(f,@sprintf("%e\n",bestIndiv[i]));
            end
            for i=1:length(searchIdx[2])
                write(f,@sprintf("%e\n",bestIndiv[i+length(searchIdx[1])]));
            end
            close(f);

            open("./FitParam/generation.dat", "w") do f
                write(f,@sprintf("%d",i));
            end
        end
        bestFitness = population[1,end];

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