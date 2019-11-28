# Residual Sum of Squares
function compute_objval_abs(
    simData::Vector{Float64},
    expData::Vector{Float64}
    )::Float64
    error::Float64 = 0.0;

    for i in eachindex(expData)
        error += (simData[i]-expData[i])^2;
    end

    return error
end


# Cosine similarity
function compute_objval_cs(
    simData::Vector{Float64},
    expData::Vector{Float64}
    )::Float64

    error::Float64 = 1.0-dot(simData,expData)/(norm(simData)*norm(expData));

    return error
end


# Define an objective function to be minimized.
function objective(
    Individual_gene::Vector{Float64},
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    SearchRegion::Matrix{Float64}
    )::Float64

    p,u0 = updateParam(
        Individual_gene::Vector{Float64},
        searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
        SearchRegion::Matrix{Float64}
    )

    # constraints --------------------------------------------------------------
    p[C.V6] = p[C.V5];
    p[C.Km6] = p[C.Km5];
    p[C.KimpDUSP] = p[C.KimDUSP];
    p[C.KexpDUSP] = p[C.KexDUSP];
    p[C.KimpcFOS] = p[C.KimFOS];
    p[C.KexpcFOS] = p[C.KexFOS];
    p[C.p52] = p[C.p47];
    p[C.m52] = p[C.m47];
    p[C.p53] = p[C.p48];
    p[C.p54] = p[C.p49];
    p[C.m54] = p[C.m49];
    p[C.p55] = p[C.p50];
    p[C.p56] = p[C.p51];
    p[C.m56] = p[C.m51];
    # --------------------------------------------------------------------------

    if Sim.simulate!(p,u0) isa Nothing
        error::Vector{Float64} = zeros(numObservables);
        for (i,target) in enumerate(observableNames)
            exp_t::Vector{Float64} = Exp.getTimepoint(i);
            normMax::Float64 = maximum(Sim.simulations[species[target],:,:]);
            if isassigned(Exp.experiments,species[target])
                error[i] = compute_objval_abs(
                    diff_sim_and_exp(
                        Sim.simulations[species[target],:,:],
                        Exp.experiments[species[target]],
                        exp_t,
                        Sim.condition,
                        normMax_sim=normMax,
                        normMax_exp=1.0
                    )...
                )
            end
        end
    
        return sum(error)
    else
        return Inf
    end
end