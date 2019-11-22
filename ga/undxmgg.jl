# Minimal Generation Gap selection for UNDX
function mggVariant(
    population::Matrix{Float64},
    n_population::Int64,
    n_children::Int64,
    n_gene::Int64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    SearchRegion::Matrix{Float64}
    )::Matrix{Float64}

    idx::Vector{Int64} = randperm(n_population);
    ip::Vector{Int64} = zeros(3);
    ip[1] = idx[1];
    ip[2] = idx[2];

    children::Matrix{Float64} = zeros(n_children,n_gene+1);
    for i = 1:n_children
        ip[3] = rand(idx[3:end]);
        children[i,:] = getNewChild(population[ip,:],n_gene,searchIdx,SearchRegion);
    end

    family::Matrix{Float64} = zeros(n_children+2,n_gene+1);
    @inbounds for i = 1:n_gene+1
        @simd for j = 1:n_children
            family[j,i] = children[j,i];
        end
        family[n_children+1,i] = population[ip[1],i];
        family[n_children+2,i] = population[ip[2],i];
    end

    family = sortslices(family,dims=1,by=x->x[end]);
    ic2::Int64 = rankSelection(n_children+2);
    @inbounds for i = 1:n_gene+1
        population[ip[1],i] = family[1,i];  # Elite
        population[ip[2],i] = family[ic2,i];  # Rank-based Roulette Selection
    end

    population = sortslices(population,dims=1,by=x->x[end]);

    return population
end


function getNewChild(
    parents::Matrix{Float64},
    n_gene::Int64,
    searchIdx::Tuple{Array{Int64,1},Array{Int64,1}},
    SearchRegion::Matrix{Float64}
    )::Vector{Float64}
    
    local child::Vector{Float64};
    MAXITER::Int8 = typemax(Int8);

    in_range::Bool = false;
    for _ in 1:MAXITER
        child = UNDX(parents,n_gene);
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

    child[end] = objective(child[1:n_gene],searchIdx,SearchRegion);

    return child
end


# Unimodal Normal Distribution Xover
function UNDX(parents::Matrix{Float64},n_gene::Int64)::Vector{Float64}
    child::Vector{Float64} = zeros(n_gene+1);
    α::Float64 = 0.5;
    β::Float64 = 0.35/sqrt(n_gene);

    p1::Vector{Float64} = [parents[1,i] for i=1:n_gene];
    p2::Vector{Float64} = [parents[2,i] for i=1:n_gene];
    p3::Vector{Float64} = [parents[3,i] for i=1:n_gene];
    d1::Float64 = norm(p2-p1);
    d2::Float64 = norm((p3-p1) - (dot((p3-p1),(p2-p1))/(d1^2))*(p2-p1));
    e1::Vector{Float64} = p1./d1;

    t::Vector{Float64} = randn(n_gene)*β*d2;
    t = t - dot(t,e1)*e1;
    t = t + randn()*α*d1*e1;

    for i = 1:n_gene
        @inbounds child[i] = t[i] + (parents[1,i]+parents[2,i])/2.0;
    end

    return child
end


function rankSelection(n_family::Int64)::Int64
    ranking::Vector{Int64} = fill(2,n_family);

    for i = 3:n_family
        ranking = append!(ranking,fill(i,n_family-i+2));
    end

    shuffle!(ranking);
    idx::Int64 = rand(1:length(ranking));

    return ranking[idx]
end