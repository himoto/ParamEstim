module ParamEstim

import .Main

using Printf
using DelimitedFiles

export
    Sim,
    search_parameter_index,
    get_search_region,
    simulate_all,
    ga_v1,
    ga_v1_continue,
    ga_v2,
    ga_v2_continue

include("ga/ga.jl")
using .GA

end  # module