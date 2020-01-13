include("ParamEstim.jl");
using .ParamEstim;

function optimize(nthParamSet::Int64)

    if !isdir("./fitparam")
        mkdir("./fitparam");
    end

    try
        files = readdir("./fitparam/$nthParamSet");
        for file in files
            if occursin(".dat",file)
                rm("./fitparam/$nthParamSet/$file");
            end
        end
    catch
        mkdir("./fitparam/$nthParamSet")
    end

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchRegion::Matrix{Float64} = getSearchRegion();

    max_generation::Int64 = 10000;
    n_population::Int64 = 5*size(searchRegion,2);
    n_children::Int64 = 50;
    n_gene::Int64 = size(searchRegion,2);
    allowable_error::Float64 = 0.35;

    (bestIndiv,bestFitness) = gaV2(
        nthParamSet,
        max_generation,
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