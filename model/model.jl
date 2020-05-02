module Model

export C, V, f_params, initial_values, diffeq

include("name2idx/name2idx.jl")
using .Name2Idx

include("set_model.jl")

end # module
