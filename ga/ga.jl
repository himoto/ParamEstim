module GA

export
    C,
    V,
    observableNames,
    numObservables,
    species,
    diff_sim_and_exp,
    Sim,
    Exp,
    f_params,
    initialValues,
    searchParameterIndex,
    getSearchRegion,
    simulateAll,
    gaV1,
    gaV1_continue,
    gaV2,
    gaV2_continue

using Printf;
using LinearAlgebra;
using Random;
using StatsBase;
using Statistics;
using DelimitedFiles;
using PyPlot;

import Seaborn;

include("../model/model.jl");
using .Model;

include("../observable.jl");

include("../experimentalData.jl");
using .Exp;

include("../simulation.jl");
using .Sim;

include("../fitness.jl");
include("../searchParameter.jl");
include("../plotFunc.jl");
include("../dynamics.jl");

include("initPop.jl");
include("converter.jl");
include("undxmgg.jl");
include("converging.jl");
include("localSearch.jl");
include("v1.jl");
include("v2.jl");

end # module