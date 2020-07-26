# Specify model parameters and/or initial values to optimize
function get_search_index()::Tuple{Array{Int64,1},Array{Int64,1}}
    # parameters
    search_idx_params::Vector{Int} = [
        C.V1
        C.Km1
        C.V5
        C.Km5
        C.V10
        C.Km10
        C.n10
        C.p11
        C.p12
        C.p13
        C.V14
        C.Km14
        C.V15
        C.Km15
        C.KimDUSP
        C.KexDUSP
        C.V20
        C.Km20
        C.V21
        C.Km21
        C.V24
        C.Km24
        C.V25
        C.Km25
        C.KimRSK
        C.KexRSK
        C.V27
        C.Km27
        C.V28
        C.Km28
        C.V29
        C.Km29
        C.V30
        C.Km30
        C.V31
        C.Km31
        C.n31
        C.p32
        C.p33
        C.p34
        C.V35
        C.Km35
        C.V36
        C.Km36
        C.V37
        C.Km37
        C.KimFOS
        C.KexFOS
        C.V42
        C.Km42
        C.V43
        C.Km43
        C.V44
        C.Km44
        C.p47
        C.m47
        C.p48
        C.p49
        C.m49
        C.p50
        C.p51
        C.m51
        C.V57
        C.Km57
        C.n57
        C.p58
        C.p59
        C.p60
        C.p61
        C.KimF
        C.KexF
        C.p63
        C.KF31
        C.nF31
        C.a
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
    #=
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
    =#

    # search_rgn[:, C.param_name] = [lower_bound, upper_bound]
    # search_rgn[:, V.var_name+length(p)] = [lower_bound, upper_bound]

    search_rgn[:, C.V1] = [7.33e-2, 6.60e-01]
    search_rgn[:, C.Km1] = [1.83e+2, 8.50e+2]
    search_rgn[:, C.V5] = [6.48e-3, 7.20e+1]
    search_rgn[:, C.Km5] = [6.00e-1, 1.60e+04]
    search_rgn[:, C.V10] = [exp(-10), exp(10)]
    search_rgn[:, C.Km10] = [exp(-10), exp(10)]
    search_rgn[:, C.n10] = [1.00, 4.00]
    search_rgn[:, C.p11] = [8.30e-13, 1.44e-2]
    search_rgn[:, C.p12] = [8.00e-8, 5.17e-2]
    search_rgn[:, C.p13] = [1.38e-7, 4.84e-1]
    search_rgn[:, C.V14] = [4.77e-3, 4.77e+1]
    search_rgn[:, C.Km14] = [2.00e+2, 2.00e+6]
    search_rgn[:, C.V15] = [exp(-10), exp(10)]
    search_rgn[:, C.Km15] = [exp(-10), exp(10)]
    search_rgn[:, C.KimDUSP] = [2.20e-4, 5.50e-1]
    search_rgn[:, C.KexDUSP] = [2.60e-4, 6.50e-1]
    search_rgn[:, C.V20] = [4.77e-3, 4.77e+1]
    search_rgn[:, C.Km20] = [2.00e+2, 2.00e+6]
    search_rgn[:, C.V21] = [exp(-10), exp(10)]
    search_rgn[:, C.Km21] = [exp(-10), exp(10)]
    search_rgn[:, C.V24] = [4.77e-2, 4.77e+0]
    search_rgn[:, C.Km24] = [2.00e+3, 2.00e+5]
    search_rgn[:, C.V25] = [exp(-10), exp(10)]
    search_rgn[:, C.Km25] = [exp(-10), exp(10)]
    search_rgn[:, C.KimRSK] = [2.20e-4, 5.50e-1]
    search_rgn[:, C.KexRSK] = [2.60e-4, 6.50e-1]
    search_rgn[:, C.V27] = [exp(-10), exp(10)]
    search_rgn[:, C.Km27] = [1.00e+2, 1.00e+4]
    search_rgn[:, C.V28] = [exp(-10), exp(10)]
    search_rgn[:, C.Km28] = [exp(-10), exp(10)]
    search_rgn[:, C.V29] = [4.77e-2, 4.77e+0]
    search_rgn[:, C.Km29] = [2.93e+3, 2.93e+5]
    search_rgn[:, C.V30] = [exp(-10), exp(10)]
    search_rgn[:, C.Km30] = [exp(-10), exp(10)]
    search_rgn[:, C.V31] = [exp(-10), exp(10)]
    search_rgn[:, C.Km31] = [exp(-10), exp(10)]
    search_rgn[:, C.n31] = [1.00, 4.00]
    search_rgn[:, C.p32] = [8.30e-13, 1.44e-2]
    search_rgn[:, C.p33] = [8.00e-8, 5.17e-2]
    search_rgn[:, C.p34] = [1.38e-7, 4.84e-1]
    search_rgn[:, C.V35] = [4.77e-3, 4.77e+1]
    search_rgn[:, C.Km35] = [2.00e+2, 2.00e+6]
    search_rgn[:, C.V36] = [exp(-10), exp(10)]
    search_rgn[:, C.Km36] = [1.00e+2, 1.00e+4]
    search_rgn[:, C.V37] = [exp(-10), exp(10)]
    search_rgn[:, C.Km37] = [exp(-10), exp(10)]
    search_rgn[:, C.KimFOS] = [2.20e-4, 5.50e-1]
    search_rgn[:, C.KexFOS] = [2.60e-4, 6.50e-1]
    search_rgn[:, C.V42] = [4.77e-3, 4.77e+1]
    search_rgn[:, C.Km42] = [2.00e+2, 2.00e+6]
    search_rgn[:, C.V43] = [exp(-10), exp(10)]
    search_rgn[:, C.Km43] = [1.00e+2, 1.00e+4]
    search_rgn[:, C.V44] = [exp(-10), exp(10)]
    search_rgn[:, C.Km44] = [exp(-10), exp(10)]
    search_rgn[:, C.p47] = [1.45e-4, 1.45e+0]
    search_rgn[:, C.m47] = [6.00e-3, 6.00e+1]
    search_rgn[:, C.p48] = [2.70e-3, 2.70e+1]
    search_rgn[:, C.p49] = [5.00e-5, 5.00e-1]
    search_rgn[:, C.m49] = [5.00e-3, 5.00e+1]
    search_rgn[:, C.p50] = [3.00e-3, 3.00e+1]
    search_rgn[:, C.p51] = [exp(-10), exp(10)]
    search_rgn[:, C.m51] = [exp(-10), exp(10)]
    search_rgn[:, C.V57] = [exp(-10), exp(10)]
    search_rgn[:, C.Km57] = [exp(-10), exp(10)]
    search_rgn[:, C.n57] = [1.00, 4.00]
    search_rgn[:, C.p58] = [8.30e-13, 1.44e-2]
    search_rgn[:, C.p59] = [8.00e-8, 5.17e-2]
    search_rgn[:, C.p60] = [1.38e-7, 4.84e-1]
    search_rgn[:, C.p61] = [exp(-10), exp(10)]
    search_rgn[:, C.KimF] = [2.20e-4, 5.50e-1]
    search_rgn[:, C.KexF] = [2.60e-4, 6.50e-1]
    search_rgn[:, C.p63] = [exp(-10), exp(10)]
    search_rgn[:, C.KF31] = [exp(-10), exp(10)]
    search_rgn[:, C.nF31] = [1.00, 4.00]
    search_rgn[:, C.a] = [1.00e+2, 5.00e+2]

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

    # constraints --------------------------------------------------------------
    p[C.V6] = p[C.V5]
    p[C.Km6] = p[C.Km5]
    p[C.KimpDUSP] = p[C.KimDUSP]
    p[C.KexpDUSP] = p[C.KexDUSP]
    p[C.KimpcFOS] = p[C.KimFOS]
    p[C.KexpcFOS] = p[C.KexFOS]
    p[C.p52] = p[C.p47]
    p[C.m52] = p[C.m47]
    p[C.p53] = p[C.p48]
    p[C.p54] = p[C.p49]
    p[C.m54] = p[C.m49]
    p[C.p55] = p[C.p50]
    p[C.p56] = p[C.p51]
    p[C.m56] = p[C.m51]
    # --------------------------------------------------------------------------

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