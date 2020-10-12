module ParamEstim

const MODEL_PATH = "models/Nakakuki_Cell_2010_ODE"

export
    C,
    V,
    observables,
    observables_index,
    Sim,
    Exp,
    Viz,
    param_values,
    initial_values,
    objective,
    get_search_index,
    get_search_region,
    decode_gene2val,
    encode_val2gene,
    encode_bestIndivVal2randGene,
    update_param,
    simulate_all,
    ga_v1,
    ga_v1_continue,
    ga_v2,
    ga_v2_continue,
    MODEL_PATH

using Printf
using LinearAlgebra
using Random
using StatsBase
using Statistics
using DelimitedFiles
using PyPlot

import Seaborn

include(MODEL_PATH * "/name2idx/parameters.jl")
include(MODEL_PATH * "/name2idx/species.jl")
include(MODEL_PATH * "/set_model.jl")
include(MODEL_PATH * "/observable.jl")
include(MODEL_PATH * "/experimental_data.jl")
include(MODEL_PATH * "/simulation.jl")
include(MODEL_PATH * "/viz.jl")
include(MODEL_PATH * "/fitness.jl")
include(MODEL_PATH * "/set_search_param.jl")

include("ga/initial_population.jl")
include("ga/undxmgg.jl")
include("ga/converging.jl")
include("ga/local_search.jl")
include("ga/v1.jl")
include("ga/v2.jl")

include("dynamics/signaling_systems.jl")
include("dynamics/temporal_dynamics.jl")

end  # module