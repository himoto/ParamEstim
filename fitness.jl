function getFitness(
    Individual_gene::Vector{Float64},
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    SearchRegion::Matrix{Float64}
    )::Float64


    p,u0 = updateParam(
        Individual_gene::Vector{Float64},
        searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
        SearchRegion::Matrix{Float64}
    )

    # constraints
    p[C.V6] = p[C.V5];
    p[C.Km6] = p[C.Km5];
    p[C.KimpDUSP] = p[C.KimDUSP];
    p[C.KexpDUSP] = p[C.KexDUSP];
    p[C.KimpcFOS] = p[C.KimFOS];
    p[C.KexpcFOS] = p[C.KexFOS];
    p[C.p52] = p[C.p47];
    p[C.m52] = p[C.m47];
    p[C.p53] = p[C.p48];
    p[C.p54] = p[C.p49];
    p[C.m54] = p[C.m49];
    p[C.p55] = p[C.p50];
    p[C.p56] = p[C.p51];
    p[C.m56] = p[C.m51];

    if Sim.numericalIntegration!(p,u0) isa Nothing
        fit::Vector{Float64} = zeros(14);

        # ERK
        norm_max::Float64 = maximum(Sim.PERK_cyt);
        fit[1] = compute_objval_abs(Sim.PERK_cyt[Int.(Exp.t2.+1),1]./norm_max,Exp.egf_ERKc_av);
        fit[2] = compute_objval_abs(Sim.PERK_cyt[Int.(Exp.t2.+1),2]./norm_max,Exp.hrg_ERKc_av);

        # RSK
        norm_max = maximum(Sim.PRSK_wcl);
        fit[3] = compute_objval_abs(Sim.PRSK_wcl[Int.(Exp.t2.+1),1]./norm_max,Exp.egf_RSKw_av);
        fit[4] = compute_objval_abs(Sim.PRSK_wcl[Int.(Exp.t2.+1),2]./norm_max,Exp.hrg_RSKw_av);

        # CREB
        norm_max = maximum(Sim.PCREB_wcl);
        fit[5] = compute_objval_abs(Sim.PCREB_wcl[Int.(Exp.t3.+1),1]./norm_max,Exp.egf_CREBw_av);
        fit[6] = compute_objval_abs(Sim.PCREB_wcl[Int.(Exp.t3.+1),2]./norm_max,Exp.hrg_CREBw_av);

        # DUSPmRNA
        norm_max = maximum(Sim.DUSPmRNA);
        fit[7] = compute_objval_abs(Sim.DUSPmRNA[Int.(Exp.t5.+1),1]./norm_max,Exp.egf_DUSPmRNA_av);
        fit[8] = compute_objval_abs(Sim.DUSPmRNA[Int.(Exp.t5.+1),2]./norm_max,Exp.hrg_DUSPmRNA_av);

        # cFosmRNA
        norm_max = maximum(Sim.cFosmRNA);
        fit[9] = compute_objval_abs(Sim.cFosmRNA[Int.(Exp.t4.+1),1]./norm_max,Exp.egf_cFosmRNA_av);
        fit[10] = compute_objval_abs(Sim.cFosmRNA[Int.(Exp.t4.+1),2]./norm_max,Exp.hrg_cFosmRNA_av);

        # cFosPro
        norm_max = maximum(Sim.cFosPro);
        fit[11] = compute_objval_abs(Sim.cFosPro[Int.(Exp.t5.+1),1]./norm_max,Exp.egf_cFosPro_av);
        fit[12] = compute_objval_abs(Sim.cFosPro[Int.(Exp.t5.+1),2]./norm_max,Exp.hrg_cFosPro_av);

        # PcFos
        norm_max = maximum(Sim.PcFos);
        fit[13] = compute_objval_abs(Sim.PcFos[Int.(Exp.t2.+1),1]./norm_max,Exp.egf_PcFos_av);
        fit[14] = compute_objval_abs(Sim.PcFos[Int.(Exp.t2.+1),2]./norm_max,Exp.hrg_PcFos_av);

        return sum(fit)
    else

        return Inf
    end
end


# Residual Sum of Squares
function compute_objval_abs(
    simData::Vector{Float64},
    expData::Vector{Float64}
    )::Float64
    fit::Float64 = 0.0;

    for i=1:length(expData)
        fit += (simData[i]-expData[i])^2;
    end

    return fit
end


# Cosine similarity
function compute_objval_cs(
    simData::Vector{Float64},
    expData::Vector{Float64}
    )::Float64

    fit::Float64 = 1.0-dot(simData,expData)/(norm(simData)*norm(expData));

    return fit
end