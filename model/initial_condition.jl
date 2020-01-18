function initial_values()::Vector{Float64}
    u0::Vector{Float64} = zeros(V.len_f_vars)

    u0[V.ERKc] = 9.60e+02
    u0[V.RSKc] = 3.53e+02
    u0[V.CREBn] = 1.00e+03
    u0[V.Elk1n] = 1.51e+03

    return u0
end