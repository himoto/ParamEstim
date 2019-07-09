function runGA()
    if !isdir("./FitParam")
        mkdir("./FitParam");
    else
        files::Vector{String} = readdir("./FitParam");
        for file in files
            rm("./FitParam/$file");
        end
    end

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchRegion::Matrix{Float64} = getSearchRegion();

    n_generation::Int64 = 1000000;
    n_population::Int64 = 5*size(searchRegion,2);
    n_children::Int64 = 50;
    n_gene::Int64 = size(searchRegion,2);
    allowable_error::Float64 = 0.35;

    (bestIndiv,bestFitness) = gaV2(
        n_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error,
        searchIdx,
        searchRegion
        );
end
