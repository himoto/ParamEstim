module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V

using Sundials

const normalization = true 
#=
if true, simulation results in each observable 
are divided by their maximum values
=#
const dt = 1.0
t = collect(0.0:dt:5400.0)  # 0, 1, 2, ..., 5400 [sec.]

const conditions = ["EGF", "HRG"]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)


function solveode(
        f::Function,
        u0::Vector{Float64},
        t::Vector{Float64},
        p::Vector{Float64})
    prob = ODEProblem(f,u0,(t[1],t[end]),p)
    try
        sol = solve(
            prob,CVODE_BDF(),
            abstol=1e-9,reltol=1e-9,dtmin=1e-8,
            saveat=dt,verbose=false
        )
        return ifelse(sol.retcode == :Success, sol, nothing)
    catch
        return nothing
    end
end


function get_steady_state!(
        f::Function,
        u0::Vector{Float64},
        p::Vector{Float64},
        eps::Float64=1e-6)::Vector{Float64}
    local sol
    while true
        sol = solveode(diffeq,u0,[0.0,dt],p)
        if sol === nothing || maximum(abs.((sol.u[end] .- u0) ./ (u0 .+ eps))) < eps
            break
        else
            for (i,val) in enumerate(sol.u[end])
                @inbounds u0[i] = val
            end
        end
    end
    return sol !== nothing ? sol.u[end] : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})
    # get steady state
    p[C.Ligand] = p[C.no_ligand]
    u0 = get_steady_state!(diffeq,u0,p)
    if isempty(u0)
        return false
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
        else
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
end
end # module