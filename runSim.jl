function data2param()::Tuple{Array{Float64,1},Array{Float64,1}}
    p::Vector{Float64} = f_params();
    u0::Vector{Float64} = initialValues();

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();

    generation::Int64 = readdlm("./FitParam/generation.dat")[1,1];
    bestIndiv::Vector{Float64} = readdlm(@sprintf("./FitParam/fitParam%d.dat",generation))[:,1];

    for i=1:length(searchIdx[1])
        p[searchIdx[1][i]] = bestIndiv[i];
    end
    for i=1:length(searchIdx[2])
        u0[searchIdx[2][i]] = bestIndiv[i+length(searchIdx[1])];
    end

    return p,u0
end


function runSim()
    local p::Vector{Float64};
    local u0::Vector{Float64};

    try
        (p,u0) = data2param();
    catch
        p = f_params();
        u0 = initialValues();
    end

    if Sim.numericalIntegration!(p,u0) isa Nothing
        plotFunc_timecourse(Sim);
    else
        println("Simulation failed.");
    end
end