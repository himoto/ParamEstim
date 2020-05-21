function simulate_all(Sim::Module;viz_type::String,show_all::Bool,stdev::Bool)
    if !isdir("./figure")
        mkdir("./figure")
    end

    if !(viz_type in ["best","average","original"])
        try
            parse(Int64,viz_type)
        catch
            error(
                "Avairable viz_type are: 'best','average','original','experiment','n(=1,2,...)'"
            )
        end
    end

    n_file::Vector{Int} = []
    if viz_type != "original"
        try
            fitparam_files::Vector{String} = readdir("./fitparam")
            for file in fitparam_files
                if occursin(r"\d",file)
                    push!(n_file, parse(Int64,file))
                end
            end
        catch
            viz_type = "original"
        end
    end

    simulaitons_all::Array{Float64,4} = fill(
        NaN,(length(observables),length(n_file),length(Sim.t),length(Sim.conditions))
    )
    if length(n_file) > 0
        if length(n_file) == 1 && viz_type == "average"
            viz_type = "best"
        end
        for (i,nth_param_set) in enumerate(n_file)
            (Sim,successful) = validate(nth_param_set)
            if successful
                for j in eachindex(observables)
                    @inbounds simulaitons_all[j,i,:,:] = Sim.simulations[j,:,:]
                end
            end
        end

        best_fitness_all::Vector{Float64} = fill(Inf,length(n_file))
        for (i,nth_param_set) in enumerate(n_file)
            if isfile("./fitparam/$nth_param_set/best_fitness.dat")
                best_fitness_all[i] = readdlm(
                    "./fitparam/$nth_param_set/best_fitness.dat"
                )[1,1]
            end
        end
        best_param_set::Int = n_file[argmin(best_fitness_all)]
        write_best_fit_param(best_param_set)

        if viz_type == "best"
            Sim,_ = validate(best_param_set)
        elseif viz_type != "average" && parse(Int64,viz_type) <= length(n_file)
            Sim,_ = validate(parse(Int64,viz_type))
        elseif viz_type != "average" && parse(Int64,viz_type) > length(n_file)
            error(
                @sprintf(
                    "n (%d) must be smaller than n_fitparam (%d)",
                    parse(Int64,viz_type), length(n_file)
                )
            )
        end

        if length(n_file) > 1
            save_param_range(n_file)
        end
    else
        p::Vector{Float64} = f_params()
        u0::Vector{Float64} = initial_values()
        if Sim.simulate!(p,u0) !== nothing
            error(
                "Simulation failed."
            )
        end
    end
    plotFunc_timecourse(
        Sim,n_file,viz_type,show_all,stdev,simulaitons_all
    )
end


function load_best_param(paramset::Int)::Tuple{Array{Float64,1},Array{Float64,1}}
    if isfile("./fitparam/$paramset/generation.dat")
        best_generation::Int64 = readdlm(
            "./fitparam/$paramset/generation.dat"
        )[1,1]
        best_indiv::Vector{Float64} = readdlm(
            @sprintf(
                "./fitparam/%d/fit_param%d.dat",
                paramset, best_generation
            )
        )[:,1]

        (p,u0) = update_param(best_indiv)
    else
        p::Vector{Float64} = f_params()
        u0::Vector{Float64} = initial_values()
    end

    return p, u0
end


function validate(nth_param_set::Int64)

    (p,u0) = load_best_param(nth_param_set)

    if Sim.simulate!(p,u0) isa Nothing
        return Sim, true
    else
        print(
            "Simulation failed.\nparameter_set #$nth_param_set"
        )
        return Sim, false
    end

    return Sim
end


function write_best_fit_param(best_param_set::Int)
    p::Vector{Float64} = f_params()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = search_parameter_index()

    best_generation::Int64 = readdlm(
        "./fitparam/$best_param_set/generation.dat"
    )[1,1]
    best_indiv = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat",
            best_param_set, best_generation
        )
    )[:,1]

    for (i,j) in enumerate(search_idx[1])
        @inbounds p[j] = best_indiv[i]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds u0[j] = best_indiv[i+length(search_idx[1])]
    end

    open("best_fit_param.txt","w") do f
        write(
            f,@sprintf(
                "# param set: %d\n",best_param_set
            )
        )
        write(
            f,"\n### Param const\n"
        )
        for i=1:C.len_f_params
            write(
                f,@sprintf(
                    "p[C.%s] = %e\n", C.param_names[i],p[i]
                )
            )
        end
        write(
            f,"\n### Non-zero initial conditions\n"
        )
        for i=1:V.len_f_vars
            if u0[i] != 0.0
                write(
                    f,@sprintf(
                        "u0[V.%s] = %e\n", V.var_names[i],u0[i]
                    )
                )
            end
        end
    end
end


function save_param_range(n_file::Vector{Int})
    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = search_parameter_index()
    popt::Matrix{Float64} = zeros(length(n_file),length(search_idx[1]))
    empty_folder::Vector{Int} = []
    for (k,nth_param_set) in enumerate(n_file)
        if !isfile("./fitparam/$nth_param_set/generation.dat")
            push!(empty_folder,k)
        else
            best_generation::Int64 = readdlm(
                "./fitparam/$nth_param_set/generation.dat"
            )[1,1]
            best_indiv = readdlm(
                @sprintf(
                    "./fitparam/%d/fit_param%d.dat",
                    nth_param_set, best_generation
                )
            )[:,1]
            popt[k,:] = best_indiv[1:length(search_idx[1])]
        end
    end
    popt = popt[setdiff(1:end,empty_folder),:]
    # --------------------------------------------------------------------------
    # Seaborn.boxplot

    fig = figure(figsize=(8,length(search_idx[1])/2.5))
    # rcParams
    rc("font",family = "Arial")
    rc("font",size = 12)
    rc("axes",linewidth = 1)
    # sns.despine
    gca().spines["right"].set_visible(false)
    gca().spines["top"].set_visible(false)
    gca().yaxis.set_ticks_position("left")
    gca().xaxis.set_ticks_position("bottom")

    ax = Seaborn.boxplot(
        data=popt,
        orient="h",
        linewidth=1,
        palette="Set2"
    )

    ax.set_xlabel("Parameter value")
    ax.set_ylabel("")
    ax.set_yticklabels([C.param_names[i] for i in search_idx[1]])
    ax.set_xscale("log")

    savefig("./figure/param_range.pdf",bbox_inches="tight")
    close(fig)
end