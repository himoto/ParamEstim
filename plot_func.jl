function plotFunc_timecourse(Sim::Module,n_file::Int64,viz_type::String,
                                show_all::Bool,stdev::Bool,simulations_all::Array{Float64,4})
    if !isdir("./figure/simulation/$viz_type")
        mkdir("./figure/simulation/$viz_type")
    end

    rc("figure",figsize = (4,3))
    rc("font",family = "Arial")
    rc("font",size = 20)
    rc("axes",linewidth = 1.2)
    rc("xtick.major",width = 1.2)
    rc("ytick.major",width = 1.2)
    rc("lines",linewidth = 2)
    rc("lines",markersize = 12)

    cmap = ["mediumblue","red"]
    shape = ["D","s"]

    for (i,name) in enumerate(observables)
        if show_all
            for j=1:n_file
                for l in eachindex(Sim.conditions)
                    plot(
                        Sim.t,simulations_all[i,j,:,l]./maximum(simulations_all[i,j,:,:]),
                        color=cmap[l],alpha=0.05
                    )
                end
            end
        end
        if viz_type != "average"
            for l in eachindex(Sim.conditions)
                plot(
                    Sim.t,Sim.simulations[i,:,l]/(maximum(Sim.simulations[i,:,:])),
                    color=cmap[l]
                )
            end
        else
            normalized = Array{Float64,4}(
                undef,length(observables),n_file,length(Sim.t),length(Sim.conditions)
            )
            for j=1:n_file
                for l in eachindex(Sim.conditions)
                    normalized[i,j,:,l] = (
                        simulations_all[i,j,:,l]./maximum(simulations_all[i,j,:,:])
                    )
                end
            end
            for l in eachindex(Sim.conditions)
                plot(
                    Sim.t,[mean(filter(!isnan,normalized[i,:,k,l])) for k in eachindex(Sim.t)],
                    color=cmap[l]
                )
            end
            if stdev
                for l in eachindex(Sim.conditions)
                    ymean = [
                        mean(filter(!isnan,normalized[i,:,k,l]))
                        for k in eachindex(Sim.t)
                    ]
                    yerr = [
                        std(filter(!isnan,normalized[i,:,k,l]),corrected=false)
                        for k in eachindex(Sim.t)
                    ]
                    fill_between(
                        Sim.t,ymean-yerr,ymean+yerr,
                        lw=0,color=cmap[l],alpha=0.1
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
                            exp_t./60.,Exp.experiments[i][condition],
                            yerr=Exp.standard_error[i][condition],
                            lw=1,markerfacecolor="None",markeredgecolor=cmap[l],ecolor=cmap[l],
                            fmt=shape[l],capsize=8,clip_on=false
                        )
                        for marker in exp_data[2]
                            marker.set_clip_on(false)
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


        xlim(0,90)
        xticks([0,30,60,90])
        yticks([0,0.3,0.6,0.9,1.2])
        ylim(0,1.2)
        xlabel("Time (min)")
        ylabel(replace(name, "_" => " "))

        savefig(
            "./figure/simulation/$viz_type/$name.pdf", bbox_inches="tight"
        )
        close()
    end
end