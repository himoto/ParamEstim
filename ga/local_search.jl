function local_search!(ip::Vector{Int64}, population::Matrix{Float64},
                        n_population::Int64, n_children::Int64, n_gene::Int64)::Matrix{Float64}
    idx::BitArray{1} = trues(n_population)
    idx[ip[1]] = false

    children::Matrix{Float64} = zeros(n_children, n_gene+1)

    for i in 1:n_children
        ip[2:end] = sample(
            collect(1:n_population)[idx], n_gene+1, replace=false
        )
        children[i,:] = mutation(
            population[ip,:], n_gene
        )
    end

    family::Matrix{Float64} = zeros(n_children+1, n_gene+1)
    @inbounds for i in 1:n_gene+1
        @simd for j in 1:n_children
            family[j,i] = children[j,i]
        end
        family[n_children+1,i] = population[ip[1],i]
    end

    family = sortslices(family, dims=1, by=x->x[end])

    for i in 1:n_gene+1
        @inbounds population[ip[1],i] = family[1,i]  # Best
    end

    population = sortslices(population, dims=1, by=x->x[end])

    return population
end


function mutation(parents::Matrix{Float64}, n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = NDM(parents, n_gene)
    for i in 1:n_gene
        @inbounds child[i] = clamp(child[i], 0.0, 1.0)
    end

    child[end] = objective(
        child[1:n_gene]
    )

    return child
end


# Normal Distribution Mutation
function NDM(parents::Matrix{Float64},n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = zeros(n_gene+1)

    GAMMA::Float64 = 0.35/sqrt(n_gene)

    t2::Vector{Float64} = zeros(n_gene)
    centroid::Vector{Float64} = reshape(
        mean(parents[2:end,1:n_gene], dims=1), n_gene
    )

    for i in 1:n_gene+1
        t2 += randn()*GAMMA*(parents[i+1,1:n_gene] - centroid)
    end

    for i in 1:n_gene
        @inbounds child[i] = parents[1,i] + t2[i]
    end

    return child
end