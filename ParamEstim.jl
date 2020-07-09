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

include("model/name2idx/parameters.jl")
include("model/name2idx/species.jl")
include("model/set_model.jl")
include("model/observable.jl")
include("model/experimental_data.jl")
include("model/simulation.jl")
include("model/visualization.jl")
include("model/fitness.jl")
include("model/set_search_param.jl")

include("ga/initial_population.jl")
include("ga/undxmgg.jl")
include("ga/converging.jl")
include("ga/local_search.jl")
include("ga/v1.jl")
include("ga/v2.jl")

include("plot_func.jl")
include("dynamics.jl")

end  # module