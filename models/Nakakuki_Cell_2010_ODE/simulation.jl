module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V

using Sundials

const STEADY_STATE_EPS = 1e-6
const MAXITER = 10
const normalization = true 
#=
if true, simulation results in each observable 
are divided by their maximum values
=#

t = collect(0.0:1.0:5400.0)  # 0, 1, 2, ..., 5400 [sec.]

const conditions = ["EGF", "HRG"]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)

function solveode(f::Function,u0::Vector{Float64},t::Vector{Float64},p::Vector{Float64})
    prob = ODEProblem(f,u0,(t[1],t[end]),p)
    try
        sol = solve(
            prob,CVODE_BDF(),
            abstol=1e-9,reltol=1e-9,dtmin=1e-8,
            saveat=1.0,verbose=false
        )
        return ifelse(length(sol.t) == length(t), sol, nothing)
    catch
        return nothing
    end
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})
    # get steady state
    p[C.Ligand] = p[C.no_ligand]
    iter::Int = 0
    while iter < MAXITER
        sol = solveode(diffeq,u0,t,p)
        if sol === nothing
            return false
        end
        if all(abs.(sol.u[end] - u0) .< STEADY_STATE_EPS)
            break
        else
            u0 .= sol.u[end]
            iter += 1
        end
    end
    # add ligand
    for (i,condition) in enumerate(conditions)
        if condition == "EGF"
            p[C.Ligand] = p[C.EGF]
        elseif condition == "HRG"
            p[C.Ligand] = p[C.HRG]
        end
        sol = solveode(diffeq,u0,t,p)
        if sol === nothing
            return false
        end
        @inbounds @simd for j in eachindex(t)
            simulations[observables_index("Phosphorylated_MEKc"),j,i] = (
                sol.u[j][V.ppMEKc]
            )
            simulations[observables_index("Phosphorylated_ERKc"),j,i] = (
                sol.u[j][V.pERKc] + sol.u[j][V.ppERKc]
            )
            simulations[observables_index("Phosphorylated_RSKw"),j,i] = (
                sol.u[j][V.pRSKc] + sol.u[j][V.pRSKn]*(p[C.Vn]/p[C.Vc])
            )
            simulations[observables_index("Phosphorylated_CREBw"),j,i] = (
                sol.u[j][V.pCREBn]*(p[C.Vn]/p[C.Vc])
            )
            simulations[observables_index("dusp_mRNA"),j,i] = (
                sol.u[j][V.duspmRNAc]
            )
            simulations[observables_index("cfos_mRNA"),j,i] = (
                sol.u[j][V.cfosmRNAc]
            )
            simulations[observables_index("cFos_Protein"),j,i] = (
                (sol.u[j][V.pcFOSn] + sol.u[j][V.cFOSn])*(p[C.Vn]/p[C.Vc])
                + sol.u[j][V.cFOSc] + sol.u[j][V.pcFOSc]
            )
            simulations[observables_index("Phosphorylated_cFos"),j,i] = (
                sol.u[j][V.pcFOSn]*(p[C.Vn]/p[C.Vc]) + sol.u[j][V.pcFOSc]
            )
        end
    end
end
end # module