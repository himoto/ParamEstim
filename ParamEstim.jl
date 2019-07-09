module ParamEstim

using Printf;

export
    runGA,
    data2param,
    runSim

include("ga/ga.jl");
using .GA;

include("runGA.jl");
include("runSim.jl");

end # module