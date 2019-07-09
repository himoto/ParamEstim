module ParamEstim

using Printf;
using DelimitedFiles;

export
    runGA,
    data2param,
    runSim

include("ga/ga.jl");
using .GA;

include("runGA.jl");
include("runSim.jl");

end # module