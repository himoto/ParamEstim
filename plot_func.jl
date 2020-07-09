function plotFunc_timecourse(Sim::Module, n_file::Vector{Int}, viz_type::String,
                                show_all::Bool, stdev::Bool, simulations_all::Array{Float64,4})
    if !isdir("./figure/simulation/$viz_type")
        mkpath("./figure/simulation/$viz_type")
    end

    rc("figure",figsize = (4,3))
    rc("font",family = "Arial")
    rc("font",size = 20)
    rc("axes",linewidth = 1.5)
    rc("xtick.major",width = 1.5)
    rc("ytick.major",width = 1.5)
    rc("lines",linewidth = 1.8)
    rc("lines",markersize = 12)

    cmap = ["mediumblue","red"]
    shape = ["D","s"]

    for (i,name) in enumerate(observables)
        gca().spines["right"].set_visible(false)
        gca().spines["top"].set_visible(false)
        gca().yaxis.set_ticks_position("left")
        gca().xaxis.set_ticks_position("bottom")
        if show_all
            for j in eachindex(n_file)
                for l in eachindex(Sim.conditions)
                    plot(
                        Sim.t,
                        simulations_all[i,j,:,l] ./ ifelse(
                            Sim.normalization,maximum(simulations_all[i,j,:,:]),1.0),
                        color=cmap[l],alpha=0.05
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
                    Sim.t,[
                        mean(
                            filter(
                                !isnan,normalized[i,:,k,l]
                            )
                        ) for k in eachindex(Sim.t)
                    ],color=cmap[l]
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
                        Sim.t,y_mean-y_std,y_mean+y_std,
                        lw=0,color=cmap[l],alpha=0.1
                    )
                end
            end
        else
            for l in eachindex(Sim.conditions)
                plot(
                    Sim.t,
                    Sim.simulations[i,:,l] / ifelse(
                        Sim.normalization,maximum(Sim.simulations[i,:,:]),1.0),
                    color=cmap[l]
                )
            end
        end

        if isassigned(Exp.experiments,i)
            exp_t = Exp.get_timepoint(i)
            if isassigned(Exp.standard_error,i)
                for (l,condition) in enumerate(Sim.conditions)
                    if condition in keys(Exp.experiments[i])
                        exp_data = errorbar(
                            exp_t./60.,Exp.experiments[i][condition],
                            yerr=Exp.standard_error[i][condition],
                            lw=1,markerfacecolor="None",
                            markeredgecolor=cmap[l],ecolor=cmap[l],
                            fmt=shape[l],capsize=8,clip_on=false
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
                            exp_t./60.,Exp.experiments[i][condition],shape[l],
                            markerfacecolor="None",markeredgecolor=cmap[l],
                            clip_on=false
                        )
                    end
                end
            end
        end

        xlim(0, 90)
        xticks([0, 30, 60, 90])
        yticks([0, 0.3, 0.6, 0.9, 1.2])
        ylim(0, 1.2)
        xlabel("Time (min)")
        ylabel(replace(replace(name, "__" => "\n"), "_" => " "))

        savefig(
            "./figure/simulation/$viz_type/$name.pdf", bbox_inches="tight"
        )
        close()
    end
end