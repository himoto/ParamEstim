module Sim
include("./observable.jl")
using ..Model

using Sundials
using SteadyStateDiffEq

const tspan = (0.0,5400.0)
const t = collect(tspan[1]:1.0:tspan[end])./60.0

const conditions = ["EGF", "HRG"]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)
function simulate!(p::Vector{Float64}, u0::Vector{Float64})
    try
        # get steady state
        p[C.Ligand] = p[C.no_ligand]
        prob = ODEProblem(diffeq,u0,(0.0,Inf),p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,DynamicSS(CVODE_BDF()),
            abstol=1e-9,reltol=1e-9,dt=1.0,verbose=false
        )
        u0 = sol.u
        # add ligand
        for (i,condition) in enumerate(conditions)
            if condition == "EGF"
                p[C.Ligand] = p[C.EGF]
            elseif condition == "HRG"
                p[C.Ligand] = p[C.HRG]
            end

            prob = ODEProblem(diffeq,u0,tspan,p)
            sol = solve(
                prob,CVODE_BDF(),saveat=1.0,dtmin=1e-8,
                abstol=1e-9,reltol=1e-9,verbose=false
            )
            @inbounds @simd for j in eachindex(t)
                simulations[obs_idx("Phosphorylated_MEKc"),j,i] = (
                    sol.u[j][V.ppMEKc]
                )
                simulations[obs_idx("Phosphorylated_ERKc"),j,i] = (
                    sol.u[j][V.pERKc] + sol.u[j][V.ppERKc]
                )
                simulations[obs_idx("Phosphorylated_RSKw"),j,i] = (
                    sol.u[j][V.pRSKc] + sol.u[j][V.pRSKn]*(p[C.Vn]/p[C.Vc])
                )
                simulations[obs_idx("Phosphorylated_CREBw"),j,i] = (
                    sol.u[j][V.pCREBn]*(p[C.Vn]/p[C.Vc])
                )
                simulations[obs_idx("dusp_mRNA"),j,i] = (
                    sol.u[j][V.duspmRNAc]
                )
                simulations[obs_idx("cfos_mRNA"),j,i] = (
                    sol.u[j][V.cfosmRNAc]
                )
                simulations[obs_idx("cFos_Protein"),j,i] = (
                    (sol.u[j][V.pcFOSn] + sol.u[j][V.cFOSn])*(p[C.Vn]/p[C.Vc])
                    + sol.u[j][V.cFOSc] + sol.u[j][V.pcFOSc]
                )
                simulations[obs_idx("Phosphorylated_cFos"),j,i] = (
                    sol.u[j][V.pcFOSn]*(p[C.Vn]/p[C.Vc]) + sol.u[j][V.pcFOSc]
                )
            end
        end
    catch
        return false
    end
end
end # module