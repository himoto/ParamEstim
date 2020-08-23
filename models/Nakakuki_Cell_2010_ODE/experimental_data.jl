module Exp
include("./observable.jl")

experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))

experiments[observables_index("Phosphorylated_MEKc")] = Dict(
    "EGF" => [0.000,0.773,0.439,0.252,0.130,0.087,0.080,0.066],
    "HRG" => [0.000,0.865,1.000,0.837,0.884,0.920,0.875,0.789],
)
error_bars[observables_index("Phosphorylated_MEKc")] = Dict(
    "EGF" => [0.000,0.030,0.048,0.009,0.009,0.017,0.012,0.008] ./ sqrt(3),
    "HRG" => [0.000,0.041,0.000,0.051,0.058,0.097,0.157,0.136] ./ sqrt(3),
)

experiments[observables_index("Phosphorylated_ERKc")] = Dict(
    "EGF" => [0.000,0.867,0.799,0.494,0.313,0.266,0.200,0.194],
    "HRG" => [0.000,0.848,1.000,0.971,0.950,0.812,0.747,0.595],
)
error_bars[observables_index("Phosphorylated_ERKc")] = Dict(
    "EGF" => [0.000,0.137,0.188,0.126,0.096,0.087,0.056,0.012] ./ sqrt(3),
    "HRG" => [0.000,0.120,0.000,0.037,0.088,0.019,0.093,0.075] ./ sqrt(3),
)

experiments[observables_index("Phosphorylated_RSKw")] = Dict(
    "EGF" => [0,0.814,0.812,0.450,0.151,0.059,0.038,0.030],
    "HRG" => [0,0.953,1.000,0.844,0.935,0.868,0.779,0.558],
)
error_bars[observables_index("Phosphorylated_RSKw")] = Dict(
    "EGF" => [0,0.064,0.194,0.030,0.027,0.031,0.043,0.051] ./ sqrt(3),
    "HRG" => [0,0.230,0.118,0.058,0.041,0.076,0.090,0.077] ./ sqrt(3),
)

experiments[observables_index("Phosphorylated_cFos")] = Dict(
    "EGF" => [0,0.060,0.109,0.083,0.068,0.049,0.027,0.017],
    "HRG" => [0,0.145,0.177,0.158,0.598,1.000,0.852,0.431],
)
error_bars[observables_index("Phosphorylated_cFos")] = Dict(
    "EGF" => [0,0.003,0.021,0.013,0.016,0.007,0.003,0.002] ./ sqrt(3),
    "HRG" => [0,0.010,0.013,0.001,0.014,0.000,0.077,0.047] ./ sqrt(3),
)

experiments[observables_index("Phosphorylated_CREBw")] = Dict(
    "EGF" => [0,0.446,0.030,0.000,0.000],
    "HRG" => [0,1.000,0.668,0.460,0.340],
)
error_bars[observables_index("Phosphorylated_CREBw")] = Dict(
    "EGF" => [0,0.0,0.0,0.0,0.0] ./ sqrt(3),
    "HRG" => [0,0.0,0.0,0.0,0.0] ./ sqrt(3),
)

experiments[observables_index("cfos_mRNA")] = Dict(
    "EGF" => [0,0.181,0.476,0.518,0.174,0.026,0.000],
    "HRG" => [0,0.353,0.861,1.000,0.637,0.300,0.059],
)
error_bars[observables_index("cfos_mRNA")] = Dict(
    "EGF" => [0.017,0.004,0.044,0.004,0.023,0.007,0.008] ./ sqrt(3),
    "HRG" => [0.017,0.006,0.065,0.044,0.087,0.023,0.001] ./ sqrt(3),
)

experiments[observables_index("cFos_Protein")] = Dict(
    "EGF" => [0,0.078,0.216,0.240,0.320,0.235],
    "HRG" => [0,0.089,0.552,0.861,1.000,0.698],
)
error_bars[observables_index("cFos_Protein")] = Dict(
    "EGF" => [0,0.036,0.028,0.056,0.071,0.048] ./ sqrt(3),
    "HRG" => [0,0.021,0.042,0.063,0.000,0.047] ./ sqrt(3),
)

experiments[observables_index("dusp_mRNA")] = Dict(
    "EGF" => [0.000,0.177,0.331,0.214,0.177,0.231],
    "HRG" => [0.000,0.221,0.750,1.000,0.960,0.934],
)
error_bars[observables_index("dusp_mRNA")] = Dict(
    "EGF" => [0.033,0.060,0.061,0.032,0.068,0.050] ./ sqrt(3),
    "HRG" => [0.027,0.059,0.094,0.124,0.113,0.108] ./ sqrt(3),
)


function get_timepoint(obs_name::String)::Vector{Float64}
    if obs_name in ["Phosphorylated_MEKc", "Phosphorylated_ERKc",
                    "Phosphorylated_RSKw", "Phosphorylated_cFos"]
        return [0., 300., 600., 900., 1800., 2700., 3600., 5400.]

    elseif obs_name == "Phosphorylated_CREBw"
        return [0., 600., 1800., 3600., 5400.]
    
    elseif obs_name == "cfos_mRNA"
        return [0.,600.,1200.,1800.,2700.,3600.,5400.]
        
    elseif obs_name in ["cFos_Protein", "dusp_mRNA"]
        return [0.,900.,1800.,2700.,3600.,5400.]
    end
end
end # module