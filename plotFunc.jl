function plotFunc_timecourse(Sim::Module);

    rc("figure",figsize = (20,8));
    rc("font",family = "Arial");
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
    plot(Sim.t,Sim.PERK_cyt[:,1]./maximum(Sim.PERK_cyt),"b");
    plot(Sim.t,Sim.PERK_cyt[:,2]./maximum(Sim.PERK_cyt),"r");
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
    plot(Sim.t,Sim.PRSK_wcl[:,1]./maximum(Sim.PRSK_wcl),"b");
    plot(Sim.t,Sim.PRSK_wcl[:,2]./maximum(Sim.PRSK_wcl),"r");
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
    plot(Sim.t,Sim.PCREB_wcl[:,1]./maximum(Sim.PCREB_wcl),"b");
    plot(Sim.t,Sim.PCREB_wcl[:,2]./maximum(Sim.PCREB_wcl),"r");
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
    plot(Sim.t,Sim.DUSPmRNA[:,1]./maximum(Sim.DUSPmRNA),"b");
    plot(Sim.t,Sim.DUSPmRNA[:,2]./maximum(Sim.DUSPmRNA),"r");
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
    plot(Sim.t,Sim.cFosmRNA[:,1]./maximum(Sim.cFosmRNA),"b");
    plot(Sim.t,Sim.cFosmRNA[:,2]./maximum(Sim.cFosmRNA),"r");
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
    plot(Sim.t,Sim.cFosPro[:,1]./maximum(Sim.cFosPro),"b");
    plot(Sim.t,Sim.cFosPro[:,2]./maximum(Sim.cFosPro),"r");
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
    plot(Sim.t,Sim.PcFos[:,1]./maximum(Sim.PcFos),"b");
    plot(Sim.t,Sim.PcFos[:,2]./maximum(Sim.PcFos),"r");
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

end