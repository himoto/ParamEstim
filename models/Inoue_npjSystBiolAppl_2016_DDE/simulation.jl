module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V

using DelayDiffEq

const normalization = true
#=
if true, simulation results in each observable
are divided by their maximum values
=#

t = collect(0.0:1.0:360.0)  # 0, 1, 2, ..., 360 [min.]
const simtime = t[end]
const sstime = 1000.0  # time to reach steady state

const conditions = ["WT"]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)

function solvedde(
        diffeq::Function,u0::Vector{Float64},history::Vector{Float64},
        tspan::Tuple{Float64,Float64},p::Vector{Float64},tau::Float64)
    h(p,t) = history
    lags = [tau]
    prob = DDEProblem(diffeq,u0,h,tspan,p;constant_lags=lags)
    alg = MethodOfSteps(BS3())
    sol = solve(
        prob,alg,saveat=1.0,abstol=1e-9,reltol=1e-9,dtmin=1e-8,verbose=false
    )
    return sol
end


function get_steady_state(
        p::Vector{Float64},u0::Vector{Float64},
        sstime::Float64,tau::Float64)::Vector{Float64}
    # get steady state (t<0)
    param::Vector{Float64} = copy(p)
    param[C.term] = 1.0
    history::Vector{Float64} = u0
    tspan::Tuple{Float64,Float64} = (0.0,sstime)
    sol = solvedde(diffeq,u0,history,tspan,param,tau)
    return sol[:,end]
end


function get_time_course(
        p::Vector{Float64},u0::Vector{Float64},
        sstime::Float64,simtime::Float64,tau::Float64)
    param::Vector{Float64} = copy(p)
    param[C.term] = 0.0
    u1::Vector{Float64} = get_steady_state(p,u0,sstime,tau)
    history::Vector{Float64} = u1
    tspan::Tuple{Float64,Float64} = (0.0,simtime)
    try
        sol = solvedde(diffeq,u1,history,tspan,param,tau)
        return ifelse(length(sol.t) == length(t), sol, nothing)
    catch
        return nothing
    end
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})
    for (i,condition) in enumerate(conditions)
        if condition == "WT"
            # pass
        end
        sol = get_time_course(p,u0,sstime,simtime,p[C.delayrnae])
        if sol === nothing
            return false
        end
        @inbounds @simd for j in eachindex(t)
            simulations[observables_index("Nuclear_NFkB"),j,i] = (
                sol.u[j][V.NFKBn]
            )
        end
    end
end
end # module