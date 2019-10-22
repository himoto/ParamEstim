function visualizeResult(Sim::Module;viz_type::String,show_all::Bool,stdev::Bool)
    if !isdir("./Fig")
        mkdir("./Fig");
    end

    if !(viz_type in ["best","average","original"])
        try
            parse(Int64,viz_type);
        catch
            error("Error: viz_type ∈ {'best','average','original',int(1~n_fitparam)}");
        end
    end

    p::Vector{Float64} = f_params();
    u0::Vector{Float64} = initialValues();

    n_file::Int = 0;
    if viz_type != "original"
        try
            fitparamFiles::Vector{String} = readdir("./FitParam");
            for file in fitparamFiles
                if occursin(r"\d",file)
                    n_file += 1;
                end
            end
        catch
            viz_type = "original";
        end
    end

    PMEK_cyt_all  = zeros((n_file,length(Sim.t),Sim.condition));
    PERK_cyt_all  = zeros((n_file,length(Sim.t),Sim.condition));
    PRSK_wcl_all  = zeros((n_file,length(Sim.t),Sim.condition));
    PCREB_wcl_all = zeros((n_file,length(Sim.t),Sim.condition));
    DUSPmRNA_all  = zeros((n_file,length(Sim.t),Sim.condition));
    cFosmRNA_all  = zeros((n_file,length(Sim.t),Sim.condition));
    cFosPro_all   = zeros((n_file,length(Sim.t),Sim.condition));
    PcFos_all     = zeros((n_file,length(Sim.t),Sim.condition));

    if n_file > 0
        if n_file == 1 && viz_type == "average"
            viz_type = "best";
        end
        for i=1:n_file
            Sim = runSimulation(i,Sim,p,u0);

            PMEK_cyt_all[i,:,:] = Sim.PMEK_cyt
            PERK_cyt_all[i,:,:] = Sim.PERK_cyt
            PRSK_wcl_all[i,:,:] = Sim.PRSK_wcl
            PCREB_wcl_all[i,:,:] = Sim.PCREB_wcl
            DUSPmRNA_all[i,:,:] = Sim.DUSPmRNA
            cFosmRNA_all[i,:,:] = Sim.cFosmRNA
            cFosPro_all[i,:,:] = Sim.cFosPro
            PcFos_all[i,:,:] = Sim.PcFos
        end

        bestFitness_all::Vector{Float64} = Inf.*ones(n_file);
        for i=1:n_file
            if isfile("./FitParam/$i/BestFitness.dat")
                bestFitness_all[i] = readdlm("./FitParam/$i/BestFitness.dat")[1,1];
            else
                bestFitness_all[i] = Inf;
            end
        end
        bestParamSet::Int = argmin(bestFitness_all);
        write_bestFitParam(bestParamSet);

        if viz_type == "best"
            Sim = runSimulation(bestParamSet,Sim,p,u0);
        elseif viz_type != "average" && parse(Int64,viz_type) <= n_file
            Sim = runSimulation(parse(Int64,viz_type),Sim,p,u0);
        elseif viz_type != "average" && parse(Int64,viz_type) > n_file
            error(@sprintf("%d is larger than n_fitparam(%d)",parse(Int64,viz_type),n_file));
        end

        if n_file > 1
            saveParamRange(n_file,p,u0);
        end

    else
        if Sim.simulate!(p,u0) !== nothing
            error("Simulation failed.");
        end

    end

    plotFunc_timecourse(Sim,n_file,viz_type,show_all,stdev,
        PMEK_cyt_all,
        PERK_cyt_all,
        PRSK_wcl_all,
        PCREB_wcl_all,
        DUSPmRNA_all,
        cFosmRNA_all,
        cFosPro_all,
        PcFos_all
    )
end


function runSimulation(nthParamSet::Int64,Sim::Module,p::Vector{Float64},u0::Vector{Float64})

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();

    # get_best_param
    try
        generation::Int64 = readdlm("./FitParam/$nthParamSet/generation.dat")[1,1];
        bestIndiv::Vector{Float64} = readdlm(@sprintf("./FitParam/%d/fitParam%d.dat",nthParamSet,generation))[:,1];

        for (i,j) in enumerate(searchIdx[1])
            p[j] = bestIndiv[i];
        end
        for (i,j) in enumerate(searchIdx[2])
            u0[j] = bestIndiv[i+length(searchIdx[1])];
        end

    catch
        # pass
    end

    if Sim.simulate!(p,u0) !== nothing
        print("Simulation failed.\nparameter_set #$nthParamSet")
    end

    return Sim
end


function saveParamRange(n_file::Int64,p::Vector{Float64},u0::Vector{Float64})
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex();
    searchParamMatrix::Matrix{Float64} = zeros(n_file,length(searchIdx[1]));

    for nthParamSet=1:n_file
        local bestIndiv::Vector{Float64};
        try
            generation::Int64 = readdlm("./FitParam/$nthParamSet/generation.dat")[1,1];
            bestIndiv = readdlm(@sprintf("./FitParam/%d/fitParam%d.dat",nthParamSet,generation))[:,1];
        catch
            bestIndiv = zeros(length(searchIdx[1])+length(searchIdx[2]));
            for (i,j) in enumerate(searchIdx[1])
                @inbounds bestIndiv[i] = p[j];
            end
            for (i,j) in enumerate(searchIdx[2])
                @inbounds bestIndiv[i+length(searchIdx[1])] = u0[j];
            end
        end
        searchParamMatrix[nthParamSet,:] = bestIndiv[1:length(searchIdx[1])];
    end

    # ==========================================================================
    # Seaborn.boxplot

    fig = figure(figsize=(8,24));
    rc("font",family = "Arial");
    rc("font",size = 12);
    rc("axes",linewidth = 1);
    gca().spines["right"].set_visible(false);
    gca().spines["top"].set_visible(false);
    gca().yaxis.set_ticks_position("left");
    gca().xaxis.set_ticks_position("bottom");

    ax = Seaborn.boxplot(
        data=searchParamMatrix,
        orient="h",
        linewidth=1,
        palette="Set2"
    );

    ax.set_xlabel("Parameter value");
    ax.set_ylabel("");
    ax.set_yticklabels([C.param_names[i] for i in searchIdx[1]]);
    ax.set_xscale("log");

    savefig("./Fig/param_range.pdf",bbox_inches="tight");
    close(fig);
    # ==========================================================================
end