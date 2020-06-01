function converging!(ip::Vector{Int64}, population::Matrix{Float64},
                        n_population::Int64, n_gene::Int64)::Matrix{Float64}
    n_children::Int8 = 10
    children::Matrix{Float64} = zeros(n_children, n_gene+1)
    for i in 1:n_children
        ip[3:end] = sample(
            collect(1:n_population), n_gene, replace=false
        )
        children[i,:] = crossover(population[ip,:], n_gene)
    end

    family::Matrix{Float64} = zeros(n_children+2, n_gene+1)
    @inbounds for i in 1:n_gene+1
        @simd for j in 1:n_children
            family[j,i] = children[j,i]
        end
        family[n_children+1,i] = population[ip[1],i]
        family[n_children+2,i] = population[ip[2],i]
    end

    family = sortslices(family, dims=1, by=x->x[end])

    ic2::Int8 = rand(2:n_children+2)
    @inbounds for i in 1:n_gene+1
        population[ip[1],i] = family[1,i]  # Best
        population[ip[2],i] = family[ic2,i]  # Random
    end

    if !isfinite(population[ip[2],end])
        population[ip[2],end] = objective(
            population[ip[2],1:n_gene]
        )
    end

    population = sortslices(population, dims=1, by=x->x[end])

    return population
end


function crossover(parents::Matrix{Float64}, n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = ENDX(parents, n_gene)
    for i in 1:n_gene
        @inbounds child[i] = clamp(child[i], 0.0, 1.0)
    end

    child[end] = Inf

    return child
end


# Extended Normal Distribution Xover
function ENDX(parents::Matrix{Float64}, n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = zeros(n_gene+1)
    
    ALPHA::Float64 = sqrt((1.0-2.0*0.35^2.0)) / 2.0
    BETA::Float64 = 0.35 / sqrt(n_gene-1)

    t1::Vector{Float64} = (parents[2,1:n_gene] - parents[1,1:n_gene])./2.0
    t2::Vector{Float64} = randn()*ALPHA*(parents[2,1:n_gene] - parents[1,1:n_gene])
    t3::Vector{Float64} = zeros(n_gene)
    centroid::Vector{Float64} = reshape(
        mean(parents[3:end,1:n_gene], dims=1), n_gene
    )

    for i in 1:n_gene
        t3 += randn()*BETA*(parents[i+2,1:n_gene] - centroid)
    end

    for i in 1:n_gene
        @inbounds child[i] = t1[i] + t2[i] + t3[i]
    end

    return child
end