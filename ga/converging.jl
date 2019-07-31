function converging(
    ip::Vector{Int64},
    population::Matrix{Float64},
    n_population::Int64,
    n_gene::Int64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    searchRegion::Matrix{Float64}
    )::Tuple{Array{Int64,1},Array{Float64,2}}
    n_children::Int8 = 10;
    children::Matrix{Float64} = zeros(n_children,n_gene+1);
    for i = 1:n_children
        ip[3:end] = sample([i for i=1:n_population],n_gene,replace=false);
        children[i,:] = crossover(population[ip,:],n_gene);
    end

    family::Matrix{Float64} = zeros(n_children+2,n_gene+1);
    for i = 1:n_gene+1
        for j = 1:n_children
            family[j,i] = children[j,i];
        end
        family[n_children+1,i] = population[ip[1],i];
        family[n_children+2,i] = population[ip[2],i];
    end

    family = sortslices(family,dims=1,by=x->x[end]);

    ic2::Int8 = rand(2:n_children+2);
    for i = 1:n_gene+1
        population[ip[1],i] = family[1,i];  # Best
        population[ip[2],i] = family[ic2,i];  # Random
    end

    if isinf(population[ip[2],end])
        population[ip[2],end] = getFitness(population[ip[2],1:n_gene],searchIdx,searchRegion);
    end

    population = sortslices(population,dims=1,by=x->x[end]);

    return ip, population
end


function crossover(parents::Matrix{Float64},n_gene::Int64)::Vector{Float64}
    local child::Vector{Float64};
    MAXITER::Int8 = typemax(Int8);
    in_range::Bool = false;
    for i = 1:MAXITER
        child = ENDX(parents,n_gene);
        if 0.0 <= minimum(child[1:n_gene]) && maximum(child[1:n_gene]) <= 1.0
            in_range = true;
            break;
        end
    end

    if !(in_range)
        for i=1:n_gene
            if child[i] < 0.0
                child[i] = 0.0;
            elseif child[i] > 1.0
                child[i] = 1.0;
            end
        end
    end

    child[end] = Inf;

    return child
end


# Extended Normal Distribution Xover
function ENDX(parents::Matrix{Float64},n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = zeros(n_gene+1);
    α::Float64 = sqrt((1.0-2.0*0.35^2.0))/2.0;
    β::Float64 = 0.35/sqrt(n_gene-1);

    p1::Vector{Float64} = parents[1,1:n_gene];
    p2::Vector{Float64} = parents[2,1:n_gene];
    t1::Vector{Float64} = (p2 - p1)./2.0;
    t2::Vector{Float64} = randn()*α*(p2-p1);
    t3::Vector{Float64} = zeros(n_gene);
    centroid::Vector{Float64} = reshape(mean(parents[3:end,1:n_gene],dims=1),n_gene);

    for i = 1:n_gene
        t3 += randn()*β*(parents[i+2,1:n_gene] - centroid);
    end

    for i = 1:n_gene
        child[i] = t1[i] + t2[i] + t3[i];
    end

    return child
end