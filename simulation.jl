module Sim
using ..GA;

using Sundials;

const tspan = (0.0,5400.0);
const t = collect(tspan[1]:1.0:tspan[end])./60.0;

const condition = 2;

PMEK_cyt  = zeros(length(t),condition);
PERK_cyt  = zeros(length(t),condition);
PRSK_wcl  = zeros(length(t),condition);
PCREB_wcl = zeros(length(t),condition);
DUSPmRNA  = zeros(length(t),condition);
cFosmRNA  = zeros(length(t),condition);
cFosPro   = zeros(length(t),condition);
PcFos     = zeros(length(t),condition);

function numericalIntegration!(p::Vector{Float64},u0::Vector{Float64})

    for i=1:condition
        if i==1
            p[C.Ligand] = p[C.EGF];
        elseif i==2
            p[C.Ligand] = p[C.HRG];
        end

        prob = ODEProblem(diffeq,u0,tspan,p);

        try
            sol = solve(prob,CVODE_BDF(),saveat=1.0,dtmin=(tspan[end]-tspan[1])/1e9,abstol=1e-9,reltol=1e-9,verbose=false);

            for j=1:length(t)
                PMEK_cyt[j,i] = sol.u[j][V.ppMEKc];
                PERK_cyt[j,i] = sol.u[j][V.pERKc] + sol.u[j][V.ppERKc];
                PRSK_wcl[j,i] = sol.u[j][V.pRSKc] + sol.u[j][V.pRSKn]*(p[C.Vn]/p[C.Vc]);
                PCREB_wcl[j,i] = sol.u[j][V.pCREBn]*(p[C.Vn]/p[C.Vc]);
                DUSPmRNA[j,i] = sol.u[j][V.duspmRNAc];
                cFosmRNA[j,i] = sol.u[j][V.cfosmRNAc];
                cFosPro[j,i] = (sol.u[j][V.pcFOSn] + sol.u[j][V.cFOSn])*(p[C.Vn]/p[C.Vc]) + sol.u[j][V.cFOSc] + sol.u[j][V.pcFOSc];
                PcFos[j,i] = sol.u[j][V.pcFOSn]*(p[C.Vn]/p[C.Vc]) + sol.u[j][V.pcFOSc];
            end
        catch

            return false
        end
    end
end
end # module