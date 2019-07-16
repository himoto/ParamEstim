function plotFunc_timecourse(Sim::Module,n_file::Int64,viz_type::String,show_all::Bool,stdev::Bool,
    PMEK_cyt_all,
    PERK_cyt_all,
    PRSK_wcl_all,
    PCREB_wcl_all,
    DUSPmRNA_all,
    cFosmRNA_all,
    cFosPro_all,
    PcFos_all
    )

    rc("figure",figsize = (20,8));
    rc("font",family = "Arial");
    rc("mathtext",fontset = "custom");
    rc("mathtext",it = "Arial:italic");
    rc("font",size = 16);
    rc("axes",linewidth = 2);
    rc("lines",linewidth = 2.5);
    rc("lines",markersize = 12);

    subplots_adjust(wspace=0.5, hspace=0.5);

    subplot(2,4,1);
    plot(Sim.t,Sim.PMEK_cyt[:,1],"b");
    plot(Sim.t,Sim.PMEK_cyt[:,2],"r");
    e = errorbar(Exp.t2/60.,Exp.egf_MEKc_av,yerr=Exp.egf_MEKc_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t2/60.,Exp.hrg_MEKc_av,yerr=Exp.hrg_MEKc_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("Phosphorylated MEK\n(cytoplasm)");

    subplot(2,4,2);
    if show_all
        for i=1:n_file
            plot(Sim.t,PERK_cyt_all[i,:,1]./maximum(PERK_cyt_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,PERK_cyt_all[i,:,2]./maximum(PERK_cyt_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.PERK_cyt[:,1]./maximum(Sim.PERK_cyt),"b");
        plot(Sim.t,Sim.PERK_cyt[:,2]./maximum(Sim.PERK_cyt),"r");
    else
        PERK_cyt_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            PERK_cyt_norm[i,:,1] = PERK_cyt_all[i,:,1]./maximum(PERK_cyt_all[i,:,:]);
            PERK_cyt_norm[i,:,2] = PERK_cyt_all[i,:,2]./maximum(PERK_cyt_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(PERK_cyt_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(PERK_cyt_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(PERK_cyt_norm[:,:,1],dims=1));
            yerr_egf = [std(PERK_cyt_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(PERK_cyt_norm[:,:,2],dims=1));
            yerr_hrg = [std(PERK_cyt_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t2/60.,Exp.egf_ERKc_av,yerr=Exp.egf_ERKc_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t2/60.,Exp.hrg_ERKc_av,yerr=Exp.hrg_ERKc_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("Phosphorylated ERK\n(cytoplasm)");

    subplot(2,4,3);
    if show_all
        for i=1:n_file
            plot(Sim.t,PRSK_wcl_all[i,:,1]./maximum(PRSK_wcl_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,PRSK_wcl_all[i,:,2]./maximum(PRSK_wcl_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.PRSK_wcl[:,1]./maximum(Sim.PRSK_wcl),"b");
        plot(Sim.t,Sim.PRSK_wcl[:,2]./maximum(Sim.PRSK_wcl),"r");
    else
        PRSK_wcl_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            PRSK_wcl_norm[i,:,1] = PRSK_wcl_all[i,:,1]./maximum(PRSK_wcl_all[i,:,:]);
            PRSK_wcl_norm[i,:,2] = PRSK_wcl_all[i,:,2]./maximum(PRSK_wcl_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(PRSK_wcl_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(PRSK_wcl_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(PRSK_wcl_norm[:,:,1],dims=1));
            yerr_egf = [std(PRSK_wcl_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(PRSK_wcl_norm[:,:,2],dims=1));
            yerr_hrg = [std(PRSK_wcl_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t2/60.,Exp.egf_RSKw_av,yerr=Exp.egf_RSKw_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t2/60.,Exp.hrg_RSKw_av,yerr=Exp.hrg_RSKw_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("Phosphorylated RSK\n(whole cell)");

    subplot(2,4,4);
    if show_all
        for i=1:n_file
            plot(Sim.t,PCREB_wcl_all[i,:,1]./maximum(PCREB_wcl_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,PCREB_wcl_all[i,:,2]./maximum(PCREB_wcl_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.PCREB_wcl[:,1]./maximum(Sim.PCREB_wcl),"b");
        plot(Sim.t,Sim.PCREB_wcl[:,2]./maximum(Sim.PCREB_wcl),"r");
    else
        PCREB_wcl_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            PCREB_wcl_norm[i,:,1] = PCREB_wcl_all[i,:,1]./maximum(PCREB_wcl_all[i,:,:]);
            PCREB_wcl_norm[i,:,2] = PCREB_wcl_all[i,:,2]./maximum(PCREB_wcl_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(PCREB_wcl_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(PCREB_wcl_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(PCREB_wcl_norm[:,:,1],dims=1));
            yerr_egf = [std(PCREB_wcl_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(PCREB_wcl_norm[:,:,2],dims=1));
            yerr_hrg = [std(PCREB_wcl_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t3/60.,Exp.egf_CREBw_av,yerr=Exp.egf_CREBw_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t3/60.,Exp.hrg_CREBw_av,yerr=Exp.hrg_CREBw_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("Phosphorylated CREB\n(whole cell)");

    subplot(2,4,5);
    if show_all
        for i=1:n_file
            plot(Sim.t,DUSPmRNA_all[i,:,1]./maximum(DUSPmRNA_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,DUSPmRNA_all[i,:,2]./maximum(DUSPmRNA_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.DUSPmRNA[:,1]./maximum(Sim.DUSPmRNA),"b");
        plot(Sim.t,Sim.DUSPmRNA[:,2]./maximum(Sim.DUSPmRNA),"r");
    else
        DUSPmRNA_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            DUSPmRNA_norm[i,:,1] = DUSPmRNA_all[i,:,1]./maximum(DUSPmRNA_all[i,:,:]);
            DUSPmRNA_norm[i,:,2] = DUSPmRNA_all[i,:,2]./maximum(DUSPmRNA_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(DUSPmRNA_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(DUSPmRNA_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(DUSPmRNA_norm[:,:,1],dims=1));
            yerr_egf = [std(DUSPmRNA_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(DUSPmRNA_norm[:,:,2],dims=1));
            yerr_hrg = [std(DUSPmRNA_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t5/60.,Exp.egf_DUSPmRNA_av,yerr=Exp.egf_DUSPmRNA_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t5/60.,Exp.hrg_DUSPmRNA_av,yerr=Exp.hrg_DUSPmRNA_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel(L"$\it{dusp}$"*" mRNA\nexpression");

    subplot(2,4,6);
    if show_all
        for i=1:n_file
            plot(Sim.t,cFosmRNA_all[i,:,1]./maximum(cFosmRNA_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,cFosmRNA_all[i,:,2]./maximum(cFosmRNA_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.cFosmRNA[:,1]./maximum(Sim.cFosmRNA),"b");
        plot(Sim.t,Sim.cFosmRNA[:,2]./maximum(Sim.cFosmRNA),"r");
    else
        cFosmRNA_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            cFosmRNA_norm[i,:,1] = cFosmRNA_all[i,:,1]./maximum(cFosmRNA_all[i,:,:]);
            cFosmRNA_norm[i,:,2] = cFosmRNA_all[i,:,2]./maximum(cFosmRNA_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(cFosmRNA_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(cFosmRNA_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(cFosmRNA_norm[:,:,1],dims=1));
            yerr_egf = [std(cFosmRNA_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(cFosmRNA_norm[:,:,2],dims=1));
            yerr_hrg = [std(cFosmRNA_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t4/60.,Exp.egf_cFosmRNA_av,yerr=Exp.egf_cFosmRNA_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t4/60.,Exp.hrg_cFosmRNA_av,yerr=Exp.hrg_cFosmRNA_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel(L"$\it{c}$"*"-"*L"$\it{fos}$"*" mRNA\nexpression");

    subplot(2,4,7);
    if show_all
        for i=1:n_file
            plot(Sim.t,cFosPro_all[i,:,1]./maximum(cFosPro_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,cFosPro_all[i,:,2]./maximum(cFosPro_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.cFosPro[:,1]./maximum(Sim.cFosPro),"b");
        plot(Sim.t,Sim.cFosPro[:,2]./maximum(Sim.cFosPro),"r");
    else
        cFosPro_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            cFosPro_norm[i,:,1] = cFosPro_all[i,:,1]./maximum(cFosPro_all[i,:,:]);
            cFosPro_norm[i,:,2] = cFosPro_all[i,:,2]./maximum(cFosPro_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(cFosPro_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(cFosPro_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(cFosPro_norm[:,:,1],dims=1));
            yerr_egf = [std(cFosPro_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(cFosPro_norm[:,:,2],dims=1));
            yerr_hrg = [std(cFosPro_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t5/60.,Exp.egf_cFosPro_av,yerr=Exp.egf_cFosPro_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t5/60.,Exp.hrg_cFosPro_av,yerr=Exp.hrg_cFosPro_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("c-Fos Protein\nexpression");

    subplot(2,4,8);
    if show_all
        for i=1:n_file
            plot(Sim.t,PcFos_all[i,:,1]./maximum(PcFos_all[i,:,:]),"b",alpha=0.05);
            plot(Sim.t,PcFos_all[i,:,2]./maximum(PcFos_all[i,:,:]),"r",alpha=0.05);
        end
    end
    if viz_type != "average"
        plot(Sim.t,Sim.PcFos[:,1]./maximum(Sim.PcFos),"b");
        plot(Sim.t,Sim.PcFos[:,2]./maximum(Sim.PcFos),"r");
    else
        PcFos_norm  = zeros(n_file,length(Sim.t),Sim.condition);
        for i=1:n_file
            PcFos_norm[i,:,1] = PcFos_all[i,:,1]./maximum(PcFos_all[i,:,:]);
            PcFos_norm[i,:,2] = PcFos_all[i,:,2]./maximum(PcFos_all[i,:,:]);
        end
        plot(Sim.t,vec(mean(PcFos_norm[:,:,1],dims=1)),"b");
        plot(Sim.t,vec(mean(PcFos_norm[:,:,2],dims=1)),"r");
        if stdev
            mean_egf = vec(mean(PcFos_norm[:,:,1],dims=1));
            yerr_egf = [std(PcFos_norm[:,i,1]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                lw=0, color="b", alpha=0.1
            );
            mean_hrg = vec(mean(PcFos_norm[:,:,2],dims=1));
            yerr_hrg = [std(PcFos_norm[:,i,2]) for i=1:length(Sim.t)];
            fill_between(
                Sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                lw=0, color="r", alpha=0.1
            );
        end
    end

    e = errorbar(Exp.t2/60.,Exp.egf_PcFos_av,yerr=Exp.egf_PcFos_se,lw=1,
        markerfacecolor="None",markeredgecolor="b",fmt="bD",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    e = errorbar(Exp.t2/60.,Exp.hrg_PcFos_av,yerr=Exp.hrg_PcFos_se,lw=1,
        markerfacecolor="None",markeredgecolor="r",fmt="rs",capsize=8,clip_on=false);
    for b in e[2]
        b.set_clip_on(false);
    end
    xlim(0,90);
    xticks([0,30,60,90]);
    yticks([0,0.2,0.4,0.6,0.8,1,1.2]);
    ylim(0,1.2);
    xlabel("Time (min)");
    ylabel("Phosphorylated c-Fos\nProtein expression");


    show();
    #savefig("./Fig/sim_$viz_type");

end