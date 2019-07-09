module GA

export
    Sim,
    Exp,
    f_params,
    initialValues,
    searchParameterIndex,
    getSearchRegion,
    plotFunc_timecourse,
    gaV1,
    gaV2

using Printf;
using LinearAlgebra;
using Random;
using StatsBase;
using Statistics;
using DelimitedFiles;
using PyPlot;

include("../model/model.jl");
using .Model;

include("../experimentalData.jl");
using .Exp;

include("../simulation.jl");
using .Sim;

include("../fitness.jl");
include("../searchParameter.jl");
include("../plotFunc.jl");

include("initPop.jl");
include("transformation.jl");
include("undxmgg.jl");
include("converging.jl");
include("localSearch.jl");
include("v1.jl");
include("v2.jl");

end # module