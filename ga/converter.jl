function decode_gene2variable(individual_gene::Vector{Float64},
                                search_region::Matrix{Float64})::Vector{Float64}
    indiv::Vector{Float64} = zeros(length(individual_gene))
    
    for i in eachindex(individual_gene)
        @inbounds indiv[i] = 10.0^(individual_gene[i]*(search_region[2,i] - search_region[1,i]) + search_region[1,i])
    end

    return round.(indiv,sigdigits=7)
end