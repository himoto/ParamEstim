function get_indiv(paramset::Int)::Vector{Float64}
    best_generation::Int64 = readdlm(
        "./fitparam/$paramset/generation.dat"
    )[1,1]
    best_indiv::Vector{Float64} = readdlm(
        @sprintf(
            "./fitparam/%d/fit_param%d.dat",
            paramset, best_generation
        )
    )[:,1]

    return best_indiv
end


function load_param(paramset::Int)::Tuple{Array{Float64,1},Array{Float64,1}}
    best_indiv::Vector{Float64} = get_indiv(paramset)
    (p,u0) = update_param(best_indiv)

    return p, u0
end


function get_executable()
    n_file::Vector{Int} = []
    fitparam_files::Vector{String} = readdir("./fitparam")
    for file in fitparam_files
        if occursin(r"\d",file)
            push!(n_file, parse(Int64,file))
        end
    end
    empty_folder::Vector{Int} = []
    for (i,nth_param_set) in enumerate(n_file)
        if !isfile("./fitparam/$nth_param_set/generation.dat")
            push!(empty_folder,i)
        end
    end
    for i in sort(empty_folder,rev=true)
        deleteat!(n_file,i)
    end

    return n_file
end


function validate(nth_param_set::Int64)

    (p,u0) = load_param(nth_param_set)

    if Sim.simulate!(p,u0) isa Nothing
        return Sim, true
    else
        print(
            "Simulation failed. #$nth_param_set\n"
        )
        return Sim, false
    end

    return Sim
end


function write_best_fit_param(best_param_set::Int)
    (p,u0) = load_param(best_param_set)

    open("best_fit_param.txt","w") do f
        write(
            f,@sprintf(
                "# param set: %d\n",best_param_set
            )
        )
        write(
            f,"\n### Param const\n"
        )
        for (i,param) in enumerate(C.NAMES)
            write(
                f,@sprintf(
                    "p[C.%s] = %e\n", param, p[i]
                )
            )
        end
        write(
            f,"\n### Non-zero initial conditions\n"
        )
        for (i,specie) in enumerate(V.NAMES)
            if u0[i] != 0.0
                write(
                    f,@sprintf(
                        "u0[V.%s] = %e\n", specie, u0[i]
                    )
                )
            end
        end
    end
end


function save_param_range(n_file::Vector{Int})
    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()
    popt::Matrix{Float64} = zeros(length(n_file),length(search_idx[1]))
    for (i,nth_param_set) in enumerate(n_file)
        best_indiv = get_indiv(nth_param_set)
        popt[i,:] = best_indiv[1:length(search_idx[1])]
    end
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
    ax.set_yticklabels([C.NAMES[i] for i in search_idx[1]])
    ax.set_xscale("log")

    savefig("./figure/param_range.pdf",bbox_inches="tight")
    close(fig)
end


function simulate_all(
        Sim::Module=Sim;
        viz_type::String,
        show_all::Bool,
        stdev::Bool)
    if !isdir("./figure")
        mkdir("./figure")
    end

    if !(viz_type in ["best","average","original","experiment"])
        try
            parse(Int64,viz_type)
        catch
            error(
                "Avairable viz_type are: 'best','average','original','experiment','n(=1,2,...)'"
            )
        end
    end

    n_file::Vector{Int} = viz_type in ["original", "experiment"] ? [] : get_executable()

    simulaitons_all::Array{Float64,4} = fill(
        NaN,(length(observables),length(n_file),length(Sim.t),length(Sim.conditions))
    )
    if viz_type != "experiment"
        if length(n_file) > 0
            if length(n_file) == 1 && viz_type == "average"
                error("viz_type should be best, not $viz_type")
                # viz_type = "best"
            end
            for (i,nth_param_set) in enumerate(n_file)
                (Sim,is_successful) = validate(nth_param_set)
                if is_successful
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
            p::Vector{Float64} = param_values()
            u0::Vector{Float64} = initial_values()
            if Sim.simulate!(p,u0) !== nothing
                error(
                    "Simulation failed."
                )
            end
        end
    end
    plotFunc_timecourse(
        Sim,n_file,viz_type,show_all,stdev,simulaitons_all
    )
end