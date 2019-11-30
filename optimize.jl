include("ParamEstim.jl");
using .ParamEstim;

function optimize(nthParamSet::Int64)

    if !isdir("./out")
        mkdir("./out");
    end

    if !isdir("./FitParam")
        mkdir("./FitParam");
    end

    try
        files = readdir("./FitParam/$nthParamSet");
        for file in files
            if occursin(".dat",file)
                rm("./FitParam/$nthParamSet/$file");
            end
        end
    catch
        mkdir("./FitParam/$nthParamSet")
    end

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchRegion::Matrix{Float64} = getSearchRegion();

    n_generation::Int64 = 10000;
    n_population::Int64 = 5*size(searchRegion,2);
    n_children::Int64 = 50;
    n_gene::Int64 = size(searchRegion,2);
    allowable_error::Float64 = 0.35;

    (bestIndiv,bestFitness) = gaV2(
        nthParamSet,
        n_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error,
        searchIdx,
        searchRegion
        );
end

###
optimize(parse(Int64,ARGS[1]));