function decode_gene2variable(individual_gene::Vector{Float64},
                                search_region::Matrix{Float64})::Vector{Float64}
    indiv_var::Vector{Float64} = zeros(length(individual_gene))
    
    for i in eachindex(individual_gene)
        @inbounds indiv_var[i] = 10^(
            individual_gene[i]*(
                search_region[2,i] - search_region[1,i]
            ) + search_region[1,i]
        )
    end

    return round.(indiv_var,sigdigits=7)
end


function encode_bestindiv2randgene(j::Int64, best_indiv::Vector{Float64}, search_region::Matrix{Float64},
                                    p0_bounds::Vector{Float64})::Float64
    
    rand_gene::Float64 = (
        log10(
            best_indiv[j]*10^(
                rand()*log10(p0_bounds[2]/p0_bounds[1])
                + log10(p0_bounds[1])
            )
        ) - search_region[1,j]
    ) / (search_region[2,j] - search_region[1,j])

    return rand_gene
end