# Specify model parameters and/or initial values to optimize
function search_parameter_index()::Tuple{Array{Int64,1},Array{Int64,1}}
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
    search_idx_initvars::Vector{Int} = [
        # V.(variableName)
    ]

    return search_idx_params, search_idx_initvars
end


function get_search_region()::Matrix{Float64}
    p::Vector{Float64} = f_params()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = search_parameter_index()
    search_param::Vector{Float64} = init_search_param(search_idx, p, u0)

    search_region::Matrix{Float64} = zeros(2, length(p)+length(u0))
    #=
    # Default: 0.1 ~ 10x
    for (i,j) in enumerate(search_idx[1])
        search_region[1,j] = search_param[i]*0.1  # lower bound
        search_region[2,j] = search_param[i]*10.0  # upper bound
    end

    # Default: 0.5 ~ 2x
    for (i,j) in enumerate(search_idx[2])
        search_region[1,j+length(p)] = search_param[i+length(search_idx[1])]*0.5  # lower bound
        search_region[2,j+length(p)] = search_param[i+length(search_idx[1])]*2.0  # upper bound
    end
    =#

    # search_region[:, C.param_name] = [lower_bound, upper_bound]
    # search_region[:, V.var_name+length(p)] = [lower_bound, upper_bound]

    search_region[:, C.V1] = [7.33e-2, 6.60e-01]
    search_region[:, C.Km1] = [1.83e+2, 8.50e+2]
    search_region[:, C.V5] = [6.48e-3, 7.20e+1]
    search_region[:, C.Km5] = [6.00e-1, 1.60e+04]
    search_region[:, C.V10] = [exp(-10), exp(10)]
    search_region[:, C.Km10] = [exp(-10), exp(10)]
    search_region[:, C.n10] = [1.00, 4.00]
    search_region[:, C.p11] = [8.30e-13, 1.44e-2]
    search_region[:, C.p12] = [8.00e-8, 5.17e-2]
    search_region[:, C.p13] = [1.38e-7, 4.84e-1]
    search_region[:, C.V14] = [4.77e-3, 4.77e+1]
    search_region[:, C.Km14] = [2.00e+2, 2.00e+6]
    search_region[:, C.V15] = [exp(-10), exp(10)]
    search_region[:, C.Km15] = [exp(-10), exp(10)]
    search_region[:, C.KimDUSP] = [2.20e-4, 5.50e-1]
    search_region[:, C.KexDUSP] = [2.60e-4, 6.50e-1]
    search_region[:, C.V20] = [4.77e-3, 4.77e+1]
    search_region[:, C.Km20] = [2.00e+2, 2.00e+6]
    search_region[:, C.V21] = [exp(-10), exp(10)]
    search_region[:, C.Km21] = [exp(-10), exp(10)]
    search_region[:, C.V24] = [4.77e-2, 4.77e+0]
    search_region[:, C.Km24] = [2.00e+3, 2.00e+5]
    search_region[:, C.V25] = [exp(-10), exp(10)]
    search_region[:, C.Km25] = [exp(-10), exp(10)]
    search_region[:, C.KimRSK] = [2.20e-4, 5.50e-1]
    search_region[:, C.KexRSK] = [2.60e-4, 6.50e-1]
    search_region[:, C.V27] = [exp(-10), exp(10)]
    search_region[:, C.Km27] = [1.00e+2, 1.00e+4]
    search_region[:, C.V28] = [exp(-10), exp(10)]
    search_region[:, C.Km28] = [exp(-10), exp(10)]
    search_region[:, C.V29] = [4.77e-2, 4.77e+0]
    search_region[:, C.Km29] = [2.93e+3, 2.93e+5]
    search_region[:, C.V30] = [exp(-10), exp(10)]
    search_region[:, C.Km30] = [exp(-10), exp(10)]
    search_region[:, C.V31] = [exp(-10), exp(10)]
    search_region[:, C.Km31] = [exp(-10), exp(10)]
    search_region[:, C.n31] = [1.00, 4.00]
    search_region[:, C.p32] = [8.30e-13, 1.44e-2]
    search_region[:, C.p33] = [8.00e-8, 5.17e-2]
    search_region[:, C.p34] = [1.38e-7, 4.84e-1]
    search_region[:, C.V35] = [4.77e-3, 4.77e+1]
    search_region[:, C.Km35] = [2.00e+2, 2.00e+6]
    search_region[:, C.V36] = [exp(-10), exp(10)]
    search_region[:, C.Km36] = [1.00e+2, 1.00e+4]
    search_region[:, C.V37] = [exp(-10), exp(10)]
    search_region[:, C.Km37] = [exp(-10), exp(10)]
    search_region[:, C.KimFOS] = [2.20e-4, 5.50e-1]
    search_region[:, C.KexFOS] = [2.60e-4, 6.50e-1]
    search_region[:, C.V42] = [4.77e-3, 4.77e+1]
    search_region[:, C.Km42] = [2.00e+2, 2.00e+6]
    search_region[:, C.V43] = [exp(-10), exp(10)]
    search_region[:, C.Km43] = [1.00e+2, 1.00e+4]
    search_region[:, C.V44] = [exp(-10), exp(10)]
    search_region[:, C.Km44] = [exp(-10), exp(10)]
    search_region[:, C.p47] = [1.45e-4, 1.45e+0]
    search_region[:, C.m47] = [6.00e-3, 6.00e+1]
    search_region[:, C.p48] = [2.70e-3, 2.70e+1]
    search_region[:, C.p49] = [5.00e-5, 5.00e-1]
    search_region[:, C.m49] = [5.00e-3, 5.00e+1]
    search_region[:, C.p50] = [3.00e-3, 3.00e+1]
    search_region[:, C.p51] = [exp(-10), exp(10)]
    search_region[:, C.m51] = [exp(-10), exp(10)]
    search_region[:, C.V57] = [exp(-10), exp(10)]
    search_region[:, C.Km57] = [exp(-10), exp(10)]
    search_region[:, C.n57] = [1.00, 4.00]
    search_region[:, C.p58] = [8.30e-13, 1.44e-2]
    search_region[:, C.p59] = [8.00e-8, 5.17e-2]
    search_region[:, C.p60] = [1.38e-7, 4.84e-1]
    search_region[:, C.p61] = [exp(-10), exp(10)]
    search_region[:, C.KimF] = [2.20e-4, 5.50e-1]
    search_region[:, C.KexF] = [2.60e-4, 6.50e-1]
    search_region[:, C.p63] = [exp(-10), exp(10)]
    search_region[:, C.KF31] = [exp(-10), exp(10)]
    search_region[:, C.nF31] = [1.00, 4.00]
    search_region[:, C.a] = [1.00e+2, 5.00e+2]

    search_region = conv_lin2log!(
        search_region, search_idx, length(p), length(search_param)
    )

    return search_region
end


function init_search_param(search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                            p::Vector{Float64}, u0::Vector{Float64})::Vector{Float64}
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
        message::String = "search_param must not contain zero."
        for (_, idx) in enumerate(search_idx[1])
            if p[idx] == 0.0
                error(
                    @sprintf(
                        "`C.%s` in search_idx_params: ", C.param_names[idx]
                    ) * message
                )
            end
        end
        for (_, idx) in enumerate(search_idx[2])
            if u0[idx] == 0.0
                error(
                    @sprintf(
                        "`V.%s` in search_idx_initvars: ", V.var_names[idx]
                    ) * message
                )
            end
        end
    end

    return search_param
end


function conv_lin2log!(search_region::Matrix{Float64}, search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                        n_param_const::Int, n_search_param::Int)::Matrix{Float64}
    for i=1:size(search_region,2)
        if minimum(search_region[:,i]) < 0.0
            message = "search_region[lb,ub] must be positive.\n"
            if i <= n_param_const
                error(
                    @sprintf(
                        "`C.%s` ", C.param_names[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.var_names[i-n_param_const]
                    ) * message
                )
            end
        elseif minimum(search_region[:,i]) == 0.0 && maximum(search_region[:,i]) != 0.0
            message = "lower_bound must be larger than 0.\n"
            if i <= n_param_const
                error(
                    @sprintf(
                        "`C.%s` ", C.param_names[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.var_names[i-n_param_const]
                    ) * message
                )
            end
        elseif search_region[2,i] - search_region[1,i] < 0.0
            message = "lower_bound < upper_bound\n"
            if i <= n_param_const
                error(
                    @sprintf(
                        "`C.%s` ", C.param_names[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.var_names[i-n_param_const]
                    ) * message
                )
            end
        end
    end

    nonzero_idx::Vector{Int} = []
    for i=1:size(search_region,2)
        if search_region[:,i] != [0.0,0.0]
            push!(nonzero_idx,i)
        end
    end
    difference::Vector{Int} = collect(
        symdiff(
            Set(nonzero_idx),
            Set(
                append!(
                    search_idx[1], n_param_const .+ search_idx[2]
                )
            )
        )
    )
    if length(difference) > 0
        for (i,j) in enumerate(difference)
            if j <= n_param_const
                print(
                    @sprintf(
                        "`C.%s`\n", C.param_names[Int(j)]
                    )
                )
            else
                print(
                    @sprintf(
                        "`V.%s`\n", V.var_names[Int(j)-n_param_const]
                    )
                )
            end
        end
        error(
            "Set these search_params in both search_idx and search_region."
        )
    end

    search_region = search_region[:,nonzero_idx]

    return log10.(search_region)
end