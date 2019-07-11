function visualizeResult(Sim::Module;viz_type::String,show_all::Bool,stdev::Bool)

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
        fitparamFiles::Vector{String} = readdir("./FitParam");
        for file in fitparamFiles
            if occursin(r"\d",file)
                n_file += 1;
            end
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

        #=
        best_fitness_all = Inf.*ones(n_file);
        for i=1:n_file
            if isfile("./FitParam/$i/BestFitness.dat")
                best_fitness_all[i] = readdlm("./FitParam/$i/BestFitness.dat");
            else
                best_fitness_all[i] = Inf;

        best_paramset = argmin(best_fitness_all);
        =#

        if viz_type == "average"
            nothing
        elseif viz_type == "best"
            nothing
            # Sim = runSimulation(Int(best_paramset),Sim,p,u0);
        elseif parse(Int64,viz_type) <= n_file
            Sim = runSimulation(parse(Int64,viz_type),Sim,p,u0);
        else
            error(@sprintf("%d is larger than n_fitparam(%d)",parse(Int64,viz_type),n_file));
        end
    else
        if !(Sim.numericalIntegration!(p,u0) isa Nothing)
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

        for i=1:length(searchIdx[1])
            p[searchIdx[1][i]] = bestIndiv[i];
        end
        for i=1:length(searchIdx[2])
            u0[searchIdx[2][i]] = bestIndiv[i+length(searchIdx[1])];
        end

    catch
        # pass
    end

    if !(Sim.numericalIntegration!(p,u0) isa Nothing)
        print("Simulation failed.\nparameter_set #$nthParamSet")
    end

    return Sim
end