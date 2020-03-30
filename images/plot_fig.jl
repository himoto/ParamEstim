function get_ylabel(name::String)
    if name == "Phosphorylated_MEKc"
        return "Phosphorylated MEK\n(cytoplasm)"
    elseif name == "Phosphorylated_ERKc"
        return "Phosphorylated ERK\n(cytoplasm)"
    elseif name == "Phosphorylated_RSKw"
        return "Phosphorylated RSK\n(whole cell)"
    elseif name == "Phosphorylated_CREBw"
        return "Phosphorylated CREB\n(whole cell)"
    elseif name == "dusp_mRNA"
        return L"$\it{dusp}$"*" mRNA\nexpression"
    elseif name == "cfos_mRNA"
        return L"$\it{c}$"*"-"*L"$\it{fos}$"*" mRNA\nexpression"
    elseif name == "cFos_Protein"
        return "c-Fos Protein\nexpression"
    elseif name == "Phosphorylated_cFos"
        return "Phosphorylated c-Fos\nProtein expression"
    end
end


function plotFunc_timecourse(Sim::Module, n_file::Vector{Int}, viz_type::String, 
                                show_all::Bool, stdev::Bool, simulations_all::Array{Float64, 4})
    if !isdir("./figure/simulation/$viz_type")
        mkdir("./figure/simulation/$viz_type")
    end

    rc("figure", figsize = (20, 8))
    rc("font", family = "Arial")
    rc("mathtext", fontset = "custom")
    rc("mathtext", it = "Arial:italic")
    rc("font", size = 18)
    rc("axes", linewidth = 2)
    rc("xtick.major", width = 2)
    rc("ytick.major", width = 2)
    rc("xtick", direction = "in")
    rc("ytick", direction = "in")
    rc("lines", linewidth = 2.5)
    rc("lines", markersize = 12)

    cmap = ["mediumblue", "red"]
    shape = ["D", "s"]

    subplots_adjust(wspace=0.5, hspace=0.5)

    for (i, name) in enumerate(observables)
        subplot(2, 4, i)
        if show_all
            for j in eachindex(n_file)
                for l in eachindex(Sim.conditions)
                    plot(
                        Sim.t, simulations_all[i, j, :, l]./maximum(simulations_all[i, j, :, :]), 
                        color=cmap[l], alpha=0.05
                    )
                end
            end
        end
        if viz_type != "average"
            for l in eachindex(Sim.conditions)
                plot(
                    Sim.t, Sim.simulations[i, :, l]/(maximum(Sim.simulations[i, :, :])), 
                    color=cmap[l]
                )
            end
        else
            normalized = Array{Float64, 4}(
                undef, length(observables), length(n_file), length(Sim.t), length(Sim.conditions)
            )
            for j in eachindex(n_file)
                for l in eachindex(Sim.conditions)
                    normalized[i, j, :, l] = simulations_all[i, j, :, l]./maximum(simulations_all[i, j, :, :])
                end
            end
            for l in eachindex(Sim.conditions)
                plot(
                    Sim.t, [mean(filter(!isnan, normalized[i, :, k, l])) for k in eachindex(Sim.t)], 
                    color=cmap[l]
                )
            end
            if stdev
                for l in eachindex(Sim.conditions)
                    ymean = [
                        mean(
                            filter(
                                !isnan, normalized[i, :, k, l]
                            )
                        ) for k in eachindex(Sim.t)
                    ]
                    yerr = [
                        std(
                            filter(
                                !isnan, normalized[i, :, k, l]
                            )
                        ) for k in eachindex(Sim.t)
                    ]
                    fill_between(
                        Sim.t, ymean - yerr, ymean + yerr, 
                        lw=0, color=cmap[l], alpha=0.1
                    )
                end
            end
        end

        if isassigned(Exp.experiments, i)
            exp_t = Exp.get_timepoint(i)
            if isassigned(Exp.standard_error, i)
                for (l, condition) in enumerate(Sim.conditions)
                    if condition in keys(Exp.experiments[i])
                        exp_data = errorbar(
                            exp_t./60., Exp.experiments[i][condition], 
                            yerr=Exp.standard_error[i][condition], 
                            ecolor=cmap[l], elinewidth=1, capsize=8, 
                            markerfacecolor="None", markeredgecolor=cmap[l], 
                            fmt=shape[l], clip_on=false
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
                for (l, condition) in enumerate(Sim.conditions)
                    if condition in keys(Exp.experiments[i])
                        plot(
                            exp_t./60., Exp.experiments[i][condition], shape[l], 
                            markerfacecolor="None", markeredgecolor=cmap[l], 
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
        ylabel(get_ylabel(name))
    end

    savefig(
        "./figure/simulation/$viz_type.png", dpi=300, bbox_inches="tight"
    )
end