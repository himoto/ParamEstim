module GA

export
    C,
    V,
    observables,
    observables_index,
    Sim,
    Exp,
    f_params,
    initial_values,
    get_search_index,
    get_search_region,
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

include("../model/model.jl")
using .Model

include("../observable.jl")

include("../experimental_data.jl")
using .Exp

include("../simulation.jl")
using .Sim

include("../fitness.jl")
include("../set_search_param.jl")
include("../plot_func.jl")
include("../dynamics.jl")

include("converter.jl")
include("initial_population.jl")
include("undxmgg.jl")
include("converging.jl")
include("local_search.jl")
include("v1.jl")
include("v2.jl")

end # module