function runGA_continue()

    nthParamSet::Int64 = parse(Int64,replace(Main.current_ipynb,r"\D"=>""));

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchRegion::Matrix{Float64} = getSearchRegion();

    n_generation::Int64 = 10000;
    n_population::Int64 = 5*size(searchRegion,2);
    n_children::Int64 = 50;
    n_gene::Int64 = size(searchRegion,2);
    allowable_error::Float64 = 0.35;

    p0_bounds = [0.1,10.0];  # [lower_bound,upper_bound]

    if !isdir("../FitParam/$nthParamSet")
        mkdir("../FitParam/$nthParamSet");
        if !isfile(@sprintf("./runGA_%d.ipynb",nthParamSet+1))
            cp(
                @sprintf("./runGA_%d.ipynb",nthParamSet),
                @sprintf("./runGA_%d.ipynb",nthParamSet+1)
            )
        end
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
    else
        (bestIndiv,bestFitness) = gaV2_continue(
            nthParamSet,
            n_generation,
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
