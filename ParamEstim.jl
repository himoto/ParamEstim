module ParamEstim

import .Main;

using Printf;
using DelimitedFiles;

export
    Sim,
    optimize,
    optimize_continue,
    visualizeResult

include("ga/ga.jl");
using .GA;

include("optimize.jl");
include("optimize_continue.jl");

end # module