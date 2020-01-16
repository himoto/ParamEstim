function simulateAll(Sim::Module;viz_type::String,show_all::Bool,stdev::Bool)
    if !isdir("./figure")
        mkdir("./figure")
    end

    if !(viz_type in ["best","average","original"])
        try
            parse(Int64,viz_type)
        catch
            error("Error: viz_type ∈ {'best','average','original',int(1~n_fitparam)}")
        end
    end

    p::Vector{Float64} = f_params()
    u0::Vector{Float64} = initialValues()

    n_file::Int = 0
    if viz_type != "original"
        try
            fitparamFiles::Vector{String} = readdir("./fitparam")
            for file in fitparamFiles
                if occursin(r"\d",file)
                    n_file += 1
                end
            end
        catch
            viz_type = "original"
        end
    end

    simulaitons_all::Array{Float64,4} = fill(
        NaN,(length(observables),n_file,length(Sim.t),length(Sim.conditions))
    )
    if n_file > 0
        if n_file == 1 && viz_type == "average"
            viz_type = "best"
        end
        for i=1:n_file
            (Sim,successful) = validate(i,p,u0)
            if successful
                for j=1:length(observables)
                    @inbounds simulaitons_all[j,i,:,:] = Sim.simulations[j,:,:]
                end
            end
        end

        bestFitness_all::Vector{Float64} = Inf.*ones(n_file)
        for i=1:n_file
            if isfile("./fitparam/$i/best_fitness.dat")
                bestFitness_all[i] = readdlm("./fitparam/$i/best_fitness.dat")[1,1]
            else
                bestFitness_all[i] = Inf
            end
        end
        bestParamSet::Int = argmin(bestFitness_all)
        write_bestFitParam(bestParamSet,p,u0)

        if viz_type == "best"
            Sim,_ = validate(bestParamSet,p,u0)
        elseif viz_type != "average" && parse(Int64,viz_type) <= n_file
            Sim,_ = validate(parse(Int64,viz_type),p,u0)
        elseif viz_type != "average" && parse(Int64,viz_type) > n_file
            error(@sprintf("%d is larger than n_fitparam(%d)",parse(Int64,viz_type),n_file))
        end

        if n_file > 1
            saveParamRange(n_file,p,u0)
        end
    else
        if Sim.simulate!(p,u0) !== nothing
            error("Simulation failed.")
        end
    end
    plotFunc_timecourse(Sim,n_file,viz_type,show_all,stdev,simulaitons_all)
end


function updateParam(paramset::Int,p::Vector{Float64},u0::Vector{Float64})

    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex()

    if isfile("./fitparam/$paramset/generation.dat")
        bestGeneration::Int64 = readdlm(
            "./fitparam/$paramset/generation.dat"
        )[1,1]
        bestIndiv::Vector{Float64} = readdlm(
            @sprintf(
                "./fitparam/%d/fit_param%d.dat",paramset,bestGeneration
            )
        )[:,1]

        for (i,j) in enumerate(searchIdx[1])
            @inbounds p[j] = bestIndiv[i]
        end
        for (i,j) in enumerate(searchIdx[2])
            @inbounds u0[j] = bestIndiv[i+length(searchIdx[1])]
        end
    end

    return p, u0
end


function validate(nthParamSet::Int64,p::Vector{Float64},u0::Vector{Float64})

    (p,u0) = updateParam(nthParamSet,p,u0)

    if Sim.simulate!(p,u0) isa Nothing
        return Sim, true
    else
        print("Simulation failed.\nparameter_set #$nthParamSet")
        return Sim, false
    end

    return Sim
end


function write_bestFitParam(bestParamSet::Int,p::Vector{Float64},u0::Vector{Float64})
    (p,u0) = updateParam(bestParamSet,p,u0)
    open("bestFitParam.txt","w") do f
        write(f,@sprintf("# param set: %d\n",bestParamSet))
        write(f,"\n### Param const\n")
        for i=1:C.len_f_params
            write(f,@sprintf("p[C.%s] = %e\n",C.param_names[i],p[i]))
        end
        write(f,"\n### Non-zero initial conditions\n")
        for i=1:V.len_f_vars
            if u0[i] != 0.0
                write(f,@sprintf("u0[V.%s] = %e\n",V.var_names[i],u0[i]))
            end
        end
    end
end


function saveParamRange(n_file::Int64,p::Vector{Float64},u0::Vector{Float64})
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}} = searchParameterIndex()
    searchParamMatrix::Matrix{Float64} = zeros(n_file,length(searchIdx[1]))

    for nthParamSet=1:n_file
        local bestIndiv::Vector{Float64}
        try
            bestGeneration::Int64 = readdlm(
                "./fitparam/$nthParamSet/generation.dat
            ")[1,1]
            bestIndiv = readdlm(
                @sprintf(
                    "./fitparam/%d/fit_param%d.dat",nthParamSet,bestGeneration
                )
            )[:,1]
        catch
            bestIndiv = zeros(length(searchIdx[1])+length(searchIdx[2]))
            for (i,j) in enumerate(searchIdx[1])
                @inbounds bestIndiv[i] = p[j]
            end
            for (i,j) in enumerate(searchIdx[2])
                @inbounds bestIndiv[i+length(searchIdx[1])] = u0[j]
            end
        end
        searchParamMatrix[nthParamSet,:] = bestIndiv[1:length(searchIdx[1])]
    end

    # --------------------------------------------------------------------------
    # Seaborn.boxplot

    fig = figure(figsize=(8,24))
    rc("font",family = "Arial")
    rc("font",size = 12)
    rc("axes",linewidth = 1)
    gca().spines["right"].set_visible(false)
    gca().spines["top"].set_visible(false)
    gca().yaxis.set_ticks_position("left")
    gca().xaxis.set_ticks_position("bottom")

    ax = Seaborn.boxplot(
        data=searchParamMatrix,
        orient="h",
        linewidth=1,
        palette="Set2"
    )

    ax.set_xlabel("Parameter value")
    ax.set_ylabel("")
    ax.set_yticklabels([C.param_names[i] for i in searchIdx[1]])
    ax.set_xscale("log")

    savefig("./figure/param_range.pdf",bbox_inches="tight")
    close(fig)
end