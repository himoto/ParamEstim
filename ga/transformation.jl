function decodeGene2Variable(
    individualGene::Vector{Float64},
    searchRegion::Matrix{Float64}
    )::Vector{Float64}
    
    X::Vector{Float64} = zeros(length(individualGene));
    
    for i in eachindex(individualGene)
        @inbounds X[i] = 10.0^(individualGene[i]*(searchRegion[2,i] - searchRegion[1,i]) + searchRegion[1,i]);
    end

    return round.(X,sigdigits=7)
end


function updateParam(
    individualGene::Vector{Float64},
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64}
    )::Tuple{Array{Float64,1},Array{Float64,1}}

    p::Vector{Float64} = f_params();
    u0::Vector{Float64} = initialValues();
    
    X::Vector{Float64} = decodeGene2Variable(individualGene,searchRegion);

    for (i,j) in enumerate(searchIdx[1])
        @inbounds p[j] = X[i];
    end
    
    for (i,j) in enumerate(searchIdx[2])
        @inbounds u0[j] = X[i+length(searchIdx[1])];
    end
    
    return p,u0
end