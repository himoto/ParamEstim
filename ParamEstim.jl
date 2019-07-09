module ParamEstim

using Printf;

export runGA, runSim

include("ga/ga.jl");
using .GA;

include("runGA.jl");
include("runSim.jl");

end # module