module Exp
include("./observable.jl")

experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))

experiments[observables_index("Nuclear_NFkB")] = Dict(
    "WT" => [
        11.505
        51.52262857
        59.289
        83.68132857
        64.46
        45.027
        34.62007143
        36.81375
        53.31616
        49.223625
        38.17835714
        39.38095714
        44.80163333
        37.23936667
        30.32244
        36.295
        35.46456
        31.04
        34.05333333
        32.175
        29.0548
        30.96566
        35.023
        28.58318571
        22.54135
    ] ./ 83.68132857
)


function get_timepoint(obs_name::String)::Vector{Float64}
    if obs_name == "Nuclear_NFkB"
        [15.0*i for i in 0:24]  # 0, 15, 30, ..., 360 [min.]
    end
end
end # module