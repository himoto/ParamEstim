module ParamEstim

import .Main;

using Printf;
using DelimitedFiles;

export
    Sim,
    runGA,
    visualizeResult

include("ga/ga.jl");
using .GA;

include("runGA.jl");

end # module