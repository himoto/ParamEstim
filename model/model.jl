module Model

export C, V, f_params, initialValues, diffeq

include("name2idx/name2idx.jl");
using .Name2Idx

include("paramConst.jl");
include("initialCondition.jl");
include("differentialEquation.jl");

end # module
