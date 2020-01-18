module Model

export C, V, f_params, initial_values, diffeq

include("name2idx/name2idx.jl");
using .Name2Idx

include("param_const.jl");
include("initial_condition.jl");
include("differential_equation.jl");

end # module
