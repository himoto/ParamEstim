module ParamEstim

import .Main;

using Printf;
using DelimitedFiles;

export
    Sim,
    searchParameterIndex,
    getSearchRegion,
    visualizeResult,
    gaV1,
    gaV1_continue,
    gaV2,
    gaV2_continue

include("ga/ga.jl");
using .GA;

end  # module