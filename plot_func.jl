function plotFunc_timecourse(Sim::Module, n_file::Vector{Int}, viz_type::String,
                                show_all::Bool, stdev::Bool, simulations_all::Array{Float64,4})
    if !isdir("./figure/simulation/$viz_type")
        mkpath("./figure/simulation/$viz_type")
    end

    Vis.set_rcParams()

    for (i,name) in enumerate(observables)
        gca().spines["right"].set_visible(false)
        gca().spines["top"].set_visible(false)
        gca().yaxis.set_ticks_position("left")
        gca().xaxis.set_ticks_position("bottom")
        if viz_type != "experiment"
            if show_all
                for j in eachindex(n_file)
                    for l in eachindex(Sim.conditions)
                        plot(
                            Sim.t ./ Viz.options[i]["divided_by"],
                            simulations_all[i,j,:,l] ./ ifelse(
                                Sim.normalization,maximum(simulations_all[i,j,:,:]),1.0),
                            color=Viz.options[i]["cmap"][l],alpha=0.05
                        )
                    end
                end
            end
            if viz_type == "average"
                normalized = Array{Float64,4}(
                    undef,length(observables),length(n_file),length(Sim.t),length(Sim.conditions)
                )
                @inbounds for j in eachindex(n_file)
                    @simd for l in eachindex(Sim.conditions)
                        normalized[i,j,:,l] = (
                            simulations_all[i,j,:,l] ./ ifelse(
                                Sim.normalization,maximum(simulations_all[i,j,:,:]),1.0
                            )
                        )
                    end
                end
                norm_max::Float64 = maximum(
                    vcat(
                        [
                            [
                                mean(
                                    filter(
                                        !isnan,normalized[i,:,k,l]
                                    )
                                ) for k in eachindex(Sim.t)
                            ] for l in eachindex(Sim.conditions)
                        ]...
                    )
                )
                for j in eachindex(n_file)
                    for k in eachindex(Sim.t)
                        for l in eachindex(Sim.conditions)
                            @inbounds normalized[i,j,k,l] /= ifelse(Sim.normalization,norm_max,1.0)
                        end
                    end
                end
                for l in eachindex(Sim.conditions)
                    plot(
                        Sim.t ./ Viz.options[i]["divided_by"],[
                            mean(
                                filter(
                                    !isnan,normalized[i,:,k,l]
                                )
                            ) for k in eachindex(Sim.t)
                        ],color=Viz.options[i]["cmap"][l]
                    )
                end
                if stdev
                    for l in eachindex(Sim.conditions)
                        y_mean = [
                            mean(
                                filter(
                                    !isnan,normalized[i,:,k,l]
                                )
                            ) for k in eachindex(Sim.t)
                        ]
                        y_std = [
                            std(
                                filter(
                                    !isnan,normalized[i,:,k,l]
                                )
                            ) for k in eachindex(Sim.t)
                        ]
                        fill_between(
                            Sim.t ./ Viz.options[i]["divided_by"],
                            y_mean-y_std,y_mean+y_std,
                            lw=0,color=Viz.options[i]["cmap"][l],alpha=0.1
                        )
                    end
                end
            else
                for l in eachindex(Sim.conditions)
                    plot(
                        Sim.t ./ Viz.options[i]["divided_by"],
                        Sim.simulations[i,:,l] / ifelse(
                            Sim.normalization,maximum(Sim.simulations[i,:,:]),1.0),
                        color=Viz.options[i]["cmap"][l]
                    )
                end
            end
        end

        if isassigned(Exp.experiments,i)
            exp_t = Exp.get_timepoint(i)
            if isassigned(Exp.standard_error,i)
                for (l,condition) in enumerate(Sim.conditions)
                    if condition in keys(Exp.experiments[i])
                        exp_data = errorbar(
                            exp_t ./ Viz.options[i]["divided_by"],
                            Exp.experiments[i][condition],
                            yerr=Exp.standard_error[i][condition],
                            lw=1,markerfacecolor="None",
                            markeredgecolor=Viz.options[i]["cmap"][l],
                            ecolor=Viz.options[i]["cmap"][l],
                            fmt=Viz.options[i]["shape"][l],capsize=8,
                            clip_on=false
                        )
                        for capline in exp_data[2]
                            capline.set_clip_on(false)
                        end
                        for barlinecol in exp_data[3]
                            barlinecol.set_clip_on(false)
                        end
                    end
                end
            else
                for (l,condition) in enumerate(Sim.conditions)
                    if condition in keys(Exp.experiments[i])
                        plot(
                            exp_t ./ Viz.options[i]["divided_by"],
                            Exp.experiments[i][condition],
                            Viz.options[i]["shape"][l],markerfacecolor="None",
                            markeredgecolor=Viz.options[i]["cmap"][l],
                            clip_on=false
                        )
                    end
                end
            end
        end

        if length(Viz.options[i]["xlim"]) > 0
            xlim(Viz.options[i]["xlim"]...)
        end
        if length(Viz.options[i]["xticks"]) > 0
            xticks(Viz.options[i]["xticks"])
        end
        if Viz.options[i]["xlabel"] !== nothing
            xlabel(Viz.options[i]["xlabel"])
        end
        if length(Viz.options[i]["ylim"]) > 0
            ylim(Viz.options[i]["ylim"]...)
        end
        if length(Viz.options[i]["yticks"]) > 0
            yticks(Viz.options[i]["yticks"])
        end
        ylabel(Viz.options[i]["ylabel"])

        savefig(
            "./figure/simulation/$viz_type/$name.pdf", bbox_inches="tight"
        )
        close()
    end
end