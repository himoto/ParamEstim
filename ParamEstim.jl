module ParamEstim

import .Main;

using Printf;
using DelimitedFiles;

export
    Sim,
    runGA,
    runGA_continue,
    visualizeResult

include("ga/ga.jl");
using .GA;

include("parest1.jl");
include("parest2.jl");

end # module