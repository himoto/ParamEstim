const observableNames = [
    "Phosphorylated_MEKc"
    "Phosphorylated_ERKc"
    "Phosphorylated_RSKw"
    "Phosphorylated_CREBw"
    "dusp_mRNA"
    "cfos_mRNA"
    "cFos_Protein"
    "Phosphorylated_cFos"
];

const numObservables = length(observableNames);
const species = Dict(zip(observableNames,1:numObservables));

const conditionNames = Dict(
    1 => "EGF",
    2 => "HRG"
)

function diff_sim_and_exp(
    simMatrix::Matrix{Float64},expDict::Dict{String,Array{Float64,1}},expTimepoint::Vector{Float64},numCondition::Int;
    normMax_sim::Float64,normMax_exp::Float64
    )
    simResult::Vector{Float64} = [];
    expResult::Vector{Float64} = [];
    
    for i in 1:numCondition
        append!(simResult,simMatrix[Int.(expTimepoint.+1),i]);
        append!(expResult,expDict[conditionNames[i]])
    end

    return (simResult./normMax_sim, expResult./normMax_exp)
end