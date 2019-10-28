module Sim
using ..Model;

using Sundials;

const tspan = (0.0,5400.0);
const t = collect(tspan[1]:1.0:tspan[end])./60.0;

const condition = 2;

PMEK_cyt  = Matrix{Float64}(undef,length(t),condition);
PERK_cyt  = Matrix{Float64}(undef,length(t),condition);
PRSK_wcl  = Matrix{Float64}(undef,length(t),condition);
PCREB_wcl = Matrix{Float64}(undef,length(t),condition);
DUSPmRNA  = Matrix{Float64}(undef,length(t),condition);
cFosmRNA  = Matrix{Float64}(undef,length(t),condition);
cFosPro   = Matrix{Float64}(undef,length(t),condition);
PcFos     = Matrix{Float64}(undef,length(t),condition);

function simulate!(p::Vector{Float64},u0::Vector{Float64})

    for i=1:condition
        if i==1
            p[C.Ligand] = p[C.EGF];
        elseif i==2
            p[C.Ligand] = p[C.HRG];
        end

        prob = ODEProblem(diffeq,u0,tspan,p);

        try
            sol = solve(prob,CVODE_BDF(),saveat=1.0,dtmin=(tspan[end]-tspan[1])/1e9,abstol=1e-9,reltol=1e-9,verbose=false);

            @inbounds @simd for j in eachindex(t)
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