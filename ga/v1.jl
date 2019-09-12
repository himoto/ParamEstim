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
    print(@sprintf("Generation%d: Best Fitness = %.6e\n",1,population[1,end]));
    flush(stdout);

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    f = open("../FitParam/$nthParamSet/fitParam1.dat", "w");
    for i=1:length(searchIdx[1])
        write(f,@sprintf("%.6e\n",bestIndiv[i]));
    end
    for i=1:length(searchIdx[2])
        write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]));
    end
    close(f);

    open("../FitParam/$nthParamSet/generation.dat", "w") do f
        write(f,@sprintf("%d",1));
    end

    open("../FitParam/$nthParamSet/bestFitness.dat", "w") do f
        write(f,@sprintf("%.6e",bestFitness));
    end

    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        bestFitness = population[1,end];
        return bestIndiv,bestFitness
    end

    for i = 2:n_generation
        population = mggVariant(population,n_population,n_children,n_gene,searchIdx,searchRegion);
        print(@sprintf("Generation%d: Best Fitness = %.6e\n",i,population[1,end]));
        flush(stdout);
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);

        if population[1,end] < bestFitness
            f = open("../FitParam/$nthParamSet/fitParam$i.dat", "w");
            for i=1:length(searchIdx[1])
                write(f,@sprintf("%.6e\n",bestIndiv[i]));
            end
            for i=1:length(searchIdx[2])
                write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]));
            end
            close(f);

            open("../FitParam/$nthParamSet/generation.dat", "w") do f
                write(f,@sprintf("%d",i));
            end
        end

        bestFitness = population[1,end];

        open("../FitParam/$nthParamSet/bestFitness.dat", "w") do f
            write(f,@sprintf("%.6e",bestFitness));
        end

        if population[1,end] <= allowable_error
            bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
            bestFitness = population[1,end];
            return bestIndiv,bestFitness
        end

        open("../FitParam/$nthParamSet/count.dat", "w") do f
            write(f,@sprintf("%d",i));
        end
    end

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    return bestIndiv,bestFitness
end


function gaV1_continue(
    nthParamSet::Int64,
    n_generation::Int64,
    n_population::Int64,
    n_children::Int64,
    n_gene::Int64,
    allowable_error::Float64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64},
    p0_bounds::Vector{Float64}
    )::Tuple{Array{Float64,1},Float64}
    
    count::Int64 = readdlm("../FitParam/$nthParamSet/count.dat")[1,1];
    generation::Int64 = readdlm("../FitParam/$nthParamSet/generation.dat")[1,1];
    bestIndiv::Vector{Float64} = readdlm(@sprintf("../FitParam/%d/fitParam%d.dat",nthParamSet,generation))[:,1];
    bestFitness::Float64 = objective(
        (log10.(bestIndiv) .- searchRegion[1,:])./(searchRegion[2,:] .- searchRegion[1,:]),
        searchIdx,
        searchRegion
    )

    population = getInitialPopulation_continue(nthParamSet,n_population,n_gene,searchIdx,searchRegion,p0_bounds);

    if bestFitness < population[1,end]
        for i=1:n_gene
            population[1,i] = (log10(bestIndiv[i])-searchRegion[1,i])/(searchRegion[2,i]-searchRegion[1,i]);
        end
        population[1,end] = bestFitness;
    else
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        bestFitness = population[1,end];
        open("../FitParam/$nthParamSet/fitParam$count.dat", "w") do f
            for i=1:n_gene
                write(f,@sprintf("%.6e",bestIndiv[i]));
            end
        end
    end

    print(@sprintf("Generation%d: Best Fitness = %.6e\n",count+1,population[1,end]));
    flush(stdout);

    if population[1,end] <= allowable_error
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
        bestFitness = population[1,end];
        return bestIndiv,bestFitness
    end

    for i = 2:n_generation
        population = mggVariant(population,n_population,n_children,n_gene,searchIdx,searchRegion);
        print(@sprintf("Generation%d: Best Fitness = %.6e\n",i+count,population[1,end]));
        flush(stdout);
        bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);

        if population[1,end] < bestFitness
            f = open(@sprintf("../FitParam/%d/fitParam%d.dat",nthParamSet,i+count), "w");
            for i=1:length(searchIdx[1])
                write(f,@sprintf("%.6e\n",bestIndiv[i]));
            end
            for i=1:length(searchIdx[2])
                write(f,@sprintf("%.6e\n",bestIndiv[i+length(searchIdx[1])]));
            end
            close(f);

            open("../FitParam/$nthParamSet/generation.dat", "w") do f
                write(f,@sprintf("%d",i+count));
            end
        end

        bestFitness = population[1,end];

        open("../FitParam/$nthParamSet/bestFitness.dat", "w") do f
            write(f,@sprintf("%.6e",bestFitness));
        end

        if population[1,end] <= allowable_error
            bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
            bestFitness = population[1,end];
            return bestIndiv,bestFitness
        end

        open("../FitParam/$nthParamSet/count.dat", "w") do f
            write(f,@sprintf("%d",i+count));
        end
    end

    bestIndiv = decodeGene2Variable(population[1,1:n_gene],searchRegion);
    bestFitness = population[1,end];

    return bestIndiv,bestFitness
end