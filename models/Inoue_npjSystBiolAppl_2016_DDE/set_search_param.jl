# Specify model parameters and/or initial values to optimize
function get_search_index()::Tuple{Array{Int64,1},Array{Int64,1}}
    # parameters
    search_idx_params::Vector{Int} = [
        C.kassaikk
        C.kassaikknfkb
        C.kassanfkbikk
        C.kdisaikk
        C.kdisaikknfkb
        C.kdisanfkbikk
        C.kassa
        C.kdisa
        C.kdegbounda
        C.kdegboundaIKK
        C.kdegfreea
        C.kdegfreeaIKK
        C.kshutboundikbain
        C.kshutboundikbaout
        C.kshutfreeikbain
        C.kshutfreeikbaout
        C.k0mrnaikba
        C.kdegmrnaikba
        C.khillprodmrnaikba
        C.kpikba
        C.kprodmrnaikba

        C.kassbikk
        C.kassbikknfkb
        C.kassbnfkbikk
        C.kdisbikk
        C.kdisbikknfkb
        C.kdisbnfkbikk
        C.kassb
        C.kdegboundb
        C.kdegboundbIKK
        C.kdegfreeb
        C.kdegfreebIKK
        C.kdisb
        C.kshutboundikbbin
        C.kshutboundikbbout
        C.kshutfreeikbbin
        C.kshutfreeikbbout
        C.k0mrnaikbb
        C.kdegmrnaikbb
        #C.khillprodmrnaikbb
        C.kpikbb
        #C.kprodmrnaikbb
        C.kasseikk
        C.kasseikknfkb
        C.kassenfkbikk
        C.kdiseikk
        C.kdiseikknfkb
        C.kdisenfkbikk

        C.kasse
        C.kdegbounde
        C.kdegboundeIKK
        C.kdegfreee
        C.kdegfreeeIKK
        C.kdise
        C.kshutboundikbein
        C.kshutboundikbeout
        C.kshutfreeikbein
        C.kshutfreeikbeout
        C.k0mrnaikbe
        C.kdegmrnaikbe
        C.khillprodmrnaikbe
        C.kpikbe
        C.kprodmrnaikbe

        C.kshutnfkbin
        C.kshutnfkbout

        C.k0mrnaa20
        C.kdegmrnaa20
        C.kprodmrnaa20
        C.kpa20
        C.kdega20
        C.khillprodmrnaa20
    ]

    # initial values
    search_idx_initials::Vector{Int} = [
        # V.(variableName)
    ]

    return search_idx_params, search_idx_initials
end


function get_search_region()::Matrix{Float64}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()
    search_param::Vector{Float64} = init_search_param(search_idx, p, u0)

    search_rgn::Matrix{Float64} = zeros(2, length(p)+length(u0))

    # Default: 0.1 ~ 10x
    for (i,j) in enumerate(search_idx[1])
        search_rgn[1,j] = search_param[i]*0.1  # lower bound
        search_rgn[2,j] = search_param[i]*10.0  # upper bound
    end

    # Default: 0.5 ~ 2x
    for (i,j) in enumerate(search_idx[2])
        search_rgn[1,j+length(p)] = search_param[i+length(search_idx[1])]*0.5  # lower bound
        search_rgn[2,j+length(p)] = search_param[i+length(search_idx[1])]*2.0  # upper bound
    end


    # search_rgn[:, C.param_name] = [lower_bound, upper_bound]
    # search_rgn[:, V.var_name+length(p)] = [lower_bound, upper_bound]
    search_rgn[:, C.khillprodmrnaikba] = [1.0, 2.0]
    search_rgn[:, C.khillprodmrnaikbe] = [1.0, 2.0]
    search_rgn[:, C.khillprodmrnaa20] = [1.0, 2.0]

    search_rgn = conv_lin2log!(search_rgn, search_idx)

    return search_rgn
end


function update_param(indiv::Vector{Float64})::Tuple{Array{Float64,1},Array{Float64,1}}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()

    for (i,j) in enumerate(search_idx[1])
        @inbounds p[j] = indiv[i]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds u0[j] = indiv[i+length(search_idx[1])]
    end

    return p, u0
end


function decode_gene2val(indiv_gene::Vector{Float64})::Vector{Float64}
    search_rgn::Matrix{Float64} = get_search_region()
    indiv::Vector{Float64} = zeros(length(indiv_gene))

    for i in eachindex(indiv_gene)
        indiv[i] = 10^(
            indiv_gene[i] * (
                search_rgn[2,i] - search_rgn[1,i]
            ) + search_rgn[1,i]
        )
    end

    return round.(indiv,sigdigits=7)
end


function encode_val2gene(indiv::Vector{Float64})
    search_rgn::Matrix{Float64} = get_search_region()
    indiv_gene::Vector{Float64} = zeros(length(indiv))

    for i in eachindex(indiv)
        indiv_gene[i] = (
            log10(indiv[i]) - search_rgn[1,i]
        ) / (
            search_rgn[2,i] - search_rgn[1,i]
        )
    end

    return indiv_gene
end


function encode_bestIndivVal2randGene(
        idx::Int64,
        best_indiv::Vector{Float64},
        p0_bounds::Vector{Float64})::Float64
    search_rgn::Matrix{Float64} = get_search_region()
    rand_gene::Float64 = (
        log10(
            best_indiv[idx]*10^(
                rand() * log10(p0_bounds[2]/p0_bounds[1]) + log10(p0_bounds[1])
            )
        ) - search_rgn[1,idx]
    ) / (
        search_rgn[2,idx] - search_rgn[1,idx]
    )
    return rand_gene
end


function init_search_param(
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
        p::Vector{Float64},
        u0::Vector{Float64})::Vector{Float64}
    duplicate::Vector{String} = []
    if length(search_idx[1]) != length(unique(search_idx[1]))
        for idx in findall([count(x->x==i,search_idx[1]) 
                            for i in unique(search_idx[1])] .!= 1)
            push!(duplicate, C.NAMES[search_idx[1][idx]])
        end
        error(
            "Duplicate parameters (C.): $duplicate"
        )
    elseif length(search_idx[2]) != length(unique(search_idx[2]))
        for idx in findall([count(x->x==i,search_idx[2]) 
                            for i in unique(search_idx[2])] .!= 1)
            push!(duplicate, V.NAMES[search_idx[2][idx]])
        end
        error(
            "Duplicate species (V.): $duplicate"
        )   
    end
    search_param = zeros(
        length(search_idx[1]) + length(search_idx[2])
    )
    for (i,j) in enumerate(search_idx[1])
        @inbounds search_param[i] = p[j]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds search_param[i+length(search_idx[1])] = u0[j]
    end

    if any(x -> x == 0.0, search_param)
        msg::String = "search_param must not contain zero."
        for idx in search_idx[1]
            if p[idx] == 0.0
                error(
                    @sprintf(
                        "`C.%s` in search_idx_params: ", C.NAMES[idx]
                    ) * msg
                )
            end
        end
        for idx in search_idx[2]
            if u0[idx] == 0.0
                error(
                    @sprintf(
                        "`V.%s` in search_idx_initials: ", V.NAMES[idx]
                    ) * msg
                )
            end
        end
    end

    return search_param
end


function conv_lin2log!(
        search_rgn::Matrix{Float64},
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}})::Matrix{Float64}
    for i=1:size(search_rgn,2)
        if minimum(search_rgn[:,i]) < 0.0
            msg = "search_rgn[lower_bound,upper_bound] must be positive.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i-C.NUM]) * msg)
            end
        elseif minimum(search_rgn[:,i]) == 0.0 && maximum(search_rgn[:,i]) != 0.0
            msg = "lower_bound must be larger than 0.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i-C.NUM]) * msg)
            end
        elseif search_rgn[2,i] - search_rgn[1,i] < 0.0
            msg = "lower_bound must be smaller than upper_bound.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i-C.NUM]) * msg)
            end
        end
    end

    nonzero_idx::Vector{Int} = []
    for i=1:size(search_rgn,2)
        if search_rgn[:,i] != [0.0,0.0]
            push!(nonzero_idx,i)
        end
    end
    difference::Vector{Int} = collect(
        symdiff(
            Set(nonzero_idx),
            Set(append!(search_idx[1], C.NUM .+ search_idx[2]))
        )
    )
    if length(difference) > 0
        for idx in difference
            if j <= C.NUM
                println(@sprintf("`C.%s`", C.NAMES[Int(idx)]))
            else
                println(@sprintf("`V.%s`", V.NAMES[Int(idx)-C.NUM]))
            end
        end
        error(
            "Set these search_params in both search_idx and search_rgn."
        )
    end

    search_rgn = search_rgn[:,nonzero_idx]

    return log10.(search_rgn)
end