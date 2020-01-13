include("ParamEstim.jl");
using .ParamEstim;

function optimize_continue(nthParamSet::Int64)

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchRegion::Matrix{Float64} = getSearchRegion();

    max_generation::Int64 = 10000;
    n_population::Int64 = 5*size(searchRegion,2);
    n_children::Int64 = 50;
    n_gene::Int64 = size(searchRegion,2);
    allowable_error::Float64 = 0.35;

    p0_bounds = [0.1,10.0];  # [lower_bound,upper_bound]

    if !isdir("./fitparam/$nthParamSet")
        mkdir("./fitparam/$nthParamSet");
        
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
    else
        (bestIndiv,bestFitness) = gaV2_continue(
            nthParamSet,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            searchIdx,
            searchRegion,
            p0_bounds
            );
    end
end

###
optimize_continue(parse(Int64,ARGS[1]));