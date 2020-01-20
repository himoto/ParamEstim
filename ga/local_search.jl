function localsearch!(ip::Vector{Int64}, population::Matrix{Float64}, n_population::Int64,
                        n_children::Int64, n_gene::Int64, search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                        search_region::Matrix{Float64})::Tuple{Array{Int64,1},Array{Float64,2}}
    idx::BitArray{1} = trues(n_population)
    idx[ip[1]] = false

    children::Matrix{Float64} = zeros(n_children,n_gene+1)

    @inbounds for i=1:n_children
        ip[2:end] = sample(
            [i for i=1:n_population][idx], n_gene+1, replace=false
        )
        children[i,:] = mutation(
            population[ip,:], n_gene, search_idx, search_region
        )
    end

    family::Matrix{Float64} = zeros(n_children+1,n_gene+1)
    @inbounds for i = 1:n_gene+1
        @simd for j = 1:n_children
            family[j,i] = children[j,i]
        end
        family[n_children+1,i] = population[ip[1],i]
    end

    family = sortslices(family, dims=1, by=x->x[end])

    for i = 1:n_gene+1
        @inbounds population[ip[1],i] = family[1,i]  # Best
    end

    population = sortslices(population, dims=1, by=x->x[end])

    return ip,population
end


function mutation(parents::Matrix{Float64}, n_gene::Int64, search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                    search_region::Matrix{Float64})::Vector{Float64}
    local child::Vector{Float64}
    MAXITER::Int8 = typemax(Int8)

    in_range::Bool = false
    for _ in 1:MAXITER
        child = NDM(parents,n_gene)
        if 0.0 <= minimum(child[1:n_gene]) && maximum(child[1:n_gene]) <= 1.0
            in_range = true
            break
        end
    end
    if !(in_range)
        for i=1:n_gene
            if child[i] < 0.0
                child[i] = 0.0
            elseif child[i] > 1.0
                child[i] = 1.0
            end
        end
    end

    child[end] = objective(
        child[1:n_gene], search_idx, search_region
    )

    return child
end


# Normal Distribution Mutation
function NDM(parents::Matrix{Float64},n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = zeros(n_gene+1)
    γ::Float64 = 0.35/sqrt(n_gene)

    p1::Vector{Float64} = parents[1,1:n_gene]
    t2::Vector{Float64} = zeros(n_gene)
    centroid::Vector{Float64} = reshape(
        mean(parents[2:end,1:n_gene], dims=1), n_gene
    )

    for i=1:n_gene+1
        t2 += randn()*γ*(parents[i+1,1:n_gene] - centroid)
    end

    for i = 1:n_gene
        @inbounds child[i] = p1[i] + t2[i]
    end

    return child
end