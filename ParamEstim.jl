module ParamEstim

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
    ga_v2_continue

using Printf
using LinearAlgebra
using Random
using StatsBase
using Statistics
using DelimitedFiles
using PyPlot

import Seaborn

include("models/Nakakuki_Cell_2010_ODE/name2idx/parameters.jl")
include("models/Nakakuki_Cell_2010_ODE/name2idx/species.jl")
include("models/Nakakuki_Cell_2010_ODE/set_model.jl")
include("models/Nakakuki_Cell_2010_ODE/observable.jl")
include("models/Nakakuki_Cell_2010_ODE/experimental_data.jl")
include("models/Nakakuki_Cell_2010_ODE/simulation.jl")
include("models/Nakakuki_Cell_2010_ODE/visualization.jl")
include("models/Nakakuki_Cell_2010_ODE/fitness.jl")
include("models/Nakakuki_Cell_2010_ODE/set_search_param.jl")

include("ga/initial_population.jl")
include("ga/undxmgg.jl")
include("ga/converging.jl")
include("ga/local_search.jl")
include("ga/v1.jl")
include("ga/v2.jl")

include("plot_func.jl")
include("dynamics.jl")

end  # module