function diffeq(du,u,p,t)
    v = Dict{Int64,Float64}()

    v[1] = p[C.V1] * p[C.a] * u[V.ppMEKc] * u[V.ERKc] /  ( p[C.Km1] * (1 + u[V.pERKc] / p[C.Km2]) + u[V.ERKc] )
    v[2] = p[C.V2] * p[C.a] * u[V.ppMEKc] * u[V.pERKc] /  ( p[C.Km2] * (1 + u[V.ERKc] / p[C.Km1]) + u[V.pERKc] )
    v[3] = p[C.V3] * u[V.pERKc] /  ( p[C.Km3] * (1 + u[V.ppERKc] / p[C.Km4]) + u[V.pERKc] )
    v[4] = p[C.V4] * u[V.ppERKc] /  ( p[C.Km4]* (1 + u[V.pERKc] / p[C.Km3]) + u[V.ppERKc] )
    v[5] = p[C.V5] * u[V.pERKn] /  ( p[C.Km5] * (1 + u[V.ppERKn] / p[C.Km6]) + u[V.pERKn] )
    v[6] = p[C.V6] * u[V.ppERKn] /  ( p[C.Km6] * (1 + u[V.pERKn] / p[C.Km5]) + u[V.ppERKn] )
    v[7] = p[C.KimERK] * u[V.ERKc] - p[C.KexERK] * (p[C.Vn]/p[C.Vc]) * u[V.ERKn]
    v[8] = p[C.KimpERK] * u[V.pERKc] - p[C.KexpERK] * (p[C.Vn]/p[C.Vc]) * u[V.pERKn]
    v[9] = p[C.KimppERK] * u[V.ppERKc] - p[C.KexppERK] * (p[C.Vn]/p[C.Vc]) * u[V.ppERKn]
    v[10] = p[C.V10] * u[V.ppERKn]^p[C.n10] / ( p[C.Km10]^p[C.n10] + u[V.ppERKn]^p[C.n10] )
    v[11] = p[C.p11] * u[V.PreduspmRNAn]
    v[12] = p[C.p12] * u[V.duspmRNAc]
    v[13] = p[C.p13] * u[V.duspmRNAc]
    v[14] = p[C.V14] * u[V.ppERKc] * u[V.DUSPc] / ( p[C.Km14] + u[V.DUSPc] )
    v[15] = p[C.V15] * u[V.pDUSPc] / ( p[C.Km15] + u[V.pDUSPc] )
    v[16] = p[C.p16] * u[V.DUSPc]
    v[17] = p[C.p17] * u[V.pDUSPc]
    v[18] = p[C.KimDUSP] * u[V.DUSPc] - p[C.KexDUSP] * (p[C.Vn]/p[C.Vc]) * u[V.DUSPn]
    v[19] = p[C.KimpDUSP] * u[V.pDUSPc] - p[C.KexpDUSP] * (p[C.Vn]/p[C.Vc]) * u[V.pDUSPn]
    v[20] = p[C.V20] * u[V.ppERKn] * u[V.DUSPn] / ( p[C.Km20] + u[V.DUSPn] )
    v[21] = p[C.V21] * u[V.pDUSPn] / ( p[C.Km21] + u[V.pDUSPn] )
    v[22] = p[C.p22] * u[V.DUSPn]
    v[23] = p[C.p23] * u[V.pDUSPn]
    v[24] = p[C.V24] * u[V.ppERKc] * u[V.RSKc] / ( p[C.Km24] + u[V.RSKc] )
    v[25] = p[C.V25] * u[V.pRSKc] / ( p[C.Km25] + u[V.pRSKc] )
    v[26] = p[C.KimRSK] * u[V.pRSKc] - p[C.KexRSK] * (p[C.Vn]/p[C.Vc]) * u[V.pRSKn]
    v[27] = p[C.V27] * u[V.pRSKn] * u[V.CREBn] / ( p[C.Km27] + u[V.CREBn] )
    v[28] = p[C.V28] * u[V.pCREBn] / ( p[C.Km28] + u[V.pCREBn] )
    v[29] = p[C.V29] * u[V.ppERKn] * u[V.Elk1n] / ( p[C.Km29] + u[V.Elk1n] )
    v[30] = p[C.V30] * u[V.pElk1n] / ( p[C.Km30] + u[V.pElk1n] )
    v[31] = p[C.V31] * (u[V.pCREBn] * u[V.pElk1n])^p[C.n31] / ( p[C.Km31]^p[C.n31] + (u[V.pCREBn] * u[V.pElk1n])^p[C.n31] + (u[V.Fn] / p[C.KF31])^p[C.nF31] )
    v[32] = p[C.p32] * u[V.PrecfosmRNAn]
    v[33] = p[C.p33] * u[V.cfosmRNAc]
    v[34] = p[C.p34] * u[V.cfosmRNAc]
    v[35] = p[C.V35] * u[V.ppERKc] * u[V.cFOSc] / ( p[C.Km35] + u[V.cFOSc] )
    v[36] = p[C.V36] * u[V.pRSKc] * u[V.cFOSc] / ( p[C.Km36] + u[V.cFOSc] )
    v[37] = p[C.V37] * u[V.pcFOSc] / ( p[C.Km37] + u[V.pcFOSc] )
    v[38] = p[C.p38] * u[V.cFOSc]
    v[39] = p[C.p39] * u[V.pcFOSc]
    v[40] = p[C.KimFOS] * u[V.cFOSc] - p[C.KexFOS] * (p[C.Vn]/p[C.Vc]) * u[V.cFOSn]
    v[41] = p[C.KimpcFOS] * u[V.pcFOSc] - p[C.KexpcFOS] * (p[C.Vn]/p[C.Vc]) * u[V.pcFOSn]
    v[42] = p[C.V42] * u[V.ppERKn] * u[V.cFOSn] / ( p[C.Km42] + u[V.cFOSn] )
    v[43] = p[C.V43] * u[V.pRSKn] * u[V.cFOSn] / ( p[C.Km43] + u[V.cFOSn] )
    v[44] = p[C.V44] * u[V.pcFOSn] / ( p[C.Km44] + u[V.pcFOSn] )
    v[45] = p[C.p45] * u[V.cFOSn]
    v[46] = p[C.p46] * u[V.pcFOSn]
    v[47] = p[C.p47] * u[V.DUSPn] * u[V.ppERKn] - p[C.m47] * u[V.DUSPn_ppERKn]
    v[48] = p[C.p48] * u[V.DUSPn_ppERKn]
    v[49] = p[C.p49] * u[V.DUSPn] * u[V.pERKn] - p[C.m49] * u[V.DUSPn_pERKn]
    v[50] = p[C.p50] * u[V.DUSPn_pERKn]
    v[51] = p[C.p51] * u[V.DUSPn] * u[V.ERKn] - p[C.m51] * u[V.DUSPn_ERKn]
    v[52] = p[C.p52] * u[V.pDUSPn] * u[V.ppERKn] - p[C.m52] * u[V.pDUSPn_ppERKn]
    v[53] = p[C.p53] * u[V.pDUSPn_ppERKn]
    v[54] = p[C.p54] * u[V.pDUSPn] * u[V.pERKn] - p[C.m54] * u[V.pDUSPn_pERKn]
    v[55] = p[C.p55] * u[V.pDUSPn_pERKn]
    v[56] = p[C.p56] * u[V.pDUSPn] * u[V.ERKn] - p[C.m56] * u[V.pDUSPn_ERKn]
    v[57] = p[C.V57] * u[V.pcFOSn]^p[C.n57] / ( p[C.Km57]^p[C.n57] + u[V.pcFOSn]^p[C.n57] )
    v[58] = p[C.p58] * u[V.PreFmRNAn]
    v[59] = p[C.p59] * u[V.FmRNAc]
    v[60] = p[C.p60] * u[V.FmRNAc]
    v[61] = p[C.p61] * u[V.Fc]
    v[62] = p[C.KimF] * u[V.Fc] - p[C.KexF] * (p[C.Vn]/p[C.Vc]) * u[V.Fn]
    v[63] = p[C.p63] * u[V.Fn]

    if p[C.Ligand] == p[C.EGF]  # EGF=10nM
        if t < 300.
            du[V.ppMEKc] = 0.00258
        elseif t < 600.
            du[V.ppMEKc] = -0.00111
        elseif t < 900.
            du[V.ppMEKc] = -0.000625
        elseif t < 1200.
            du[V.ppMEKc] = -0.000135
        elseif t < 1800.
            du[V.ppMEKc] = -0.000135
        elseif t < 2700.
            du[V.ppMEKc] = -0.0000480
        elseif t < 3600.
            du[V.ppMEKc] = -0.00000852
        elseif t <= 5400.
            du[V.ppMEKc] = -0.00000728
        end
    elseif p[C.Ligand] == p[C.HRG]  # HRG=10nM
        if t < 300.
            du[V.ppMEKc] = 0.00288
        elseif t < 600.
            du[V.ppMEKc] = 0.000451
        elseif t < 900.
            du[V.ppMEKc] = -0.000545
        elseif t < 1200.
            du[V.ppMEKc] = 0.0000522
        elseif t < 1800.
            du[V.ppMEKc] = 0.0000522
        elseif t < 2700.
            du[V.ppMEKc] = 0.0000399
        elseif t < 3600.
            du[V.ppMEKc] = -0.0000500
        elseif t <= 5400.
            du[V.ppMEKc] = -0.0000478
        end
    else
        du[V.ppMEKc] = 0.0
    end

    du[V.CREBn] = -v[27] + v[28]
    du[V.pCREBn] = v[27] - v[28]
    du[V.ERKc] = -v[1] + v[3] - v[7]
    du[V.ERKn] = v[5] + v[7]*(p[C.Vc]/p[C.Vn]) + v[50] -v[51] + v[55] -v[56]
    du[V.pERKc] = v[1] - v[2] -v[3] +v[4] -v[8]
    du[V.pERKn] = -v[5] + v[6] + v[8]*(p[C.Vc]/p[C.Vn]) + v[48] - v[49] + v[53] - v[54]
    du[V.ppERKc] =  v[2] - v[4] - v[9]
    du[V.ppERKn] = -v[6] + v[9]*(p[C.Vc]/p[C.Vn]) - v[47] - v[52]
    du[V.Elk1n] = -v[29] + v[30]
    du[V.pElk1n] = v[29] - v[30]
    du[V.cFOSc] = v[34] - v[35] - v[36] + v[37] - v[38] - v[40]
    du[V.cFOSn] = v[40]*(p[C.Vc]/p[C.Vn]) - v[42] - v[43] + v[44] - v[45]
    du[V.pcFOSc] = v[35] + v[36] - v[37] - v[39] - v[41]
    du[V.pcFOSn] = v[41]*(p[C.Vc]/p[C.Vn]) + v[42] + v[43] - v[44] - v[46]
    du[V.DUSPc] = v[13] - v[14] + v[15] - v[16] - v[18]
    du[V.DUSPn] = v[18]*(p[C.Vc]/p[C.Vn]) - v[20] + v[21] - v[22] - v[47] + v[48] - v[49] + v[50] - v[51]
    du[V.pDUSPc] = v[14] - v[15] - v[17] - v[19]
    du[V.pDUSPn] = v[19]*(p[C.Vc]/p[C.Vn]) + v[20] - v[21] - v[23] - v[52] + v[53] - v[54] + v[55] - v[56]
    du[V.DUSPn_ERKn] = v[51]
    du[V.DUSPn_pERKn] = v[49] - v[50]
    du[V.DUSPn_ppERKn] = v[47] - v[48]
    du[V.pDUSPn_ERKn] = v[56]
    du[V.pDUSPn_pERKn] = v[54] - v[55]
    du[V.pDUSPn_ppERKn] = v[52] - v[53]
    du[V.RSKc] = -v[24] + v[25]
    du[V.pRSKc] = v[24] - v[25] - v[26]
    du[V.pRSKn] = v[26]*(p[C.Vc]/p[C.Vn])
    du[V.PrecfosmRNAn] = v[31] - v[32]
    du[V.PreduspmRNAn] = v[10] - v[11]
    du[V.cfosmRNAc] = v[32]*(p[C.Vn]/p[C.Vc]) - v[33]
    du[V.duspmRNAc] = v[11]*(p[C.Vn]/p[C.Vc]) - v[12]
    du[V.Fc] = v[60] - v[61] - v[62]
    du[V.Fn] = v[62]*(p[C.Vc]/p[C.Vn]) - v[63]
    du[V.FmRNAc] = v[58]*(p[C.Vn]/p[C.Vc]) - v[59]
    du[V.PreFmRNAn] = v[57] - v[58]

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = zeros(C.NUM)

    p[C.V1] = 0.34284837
    p[C.Km1] = 307.0415253
    p[C.V2] = 2.20e-01
    p[C.Km2] = 3.50e+02
    p[C.V3] = 7.20e-01
    p[C.Km3] = 1.60e+02
    p[C.V4] = 6.48e-01
    p[C.Km4] = 6.00e+01
    p[C.V5] = 19.49872346
    p[C.Km5]  =29.94073716
    p[C.V6] = p[C.V5]
    p[C.Km6] = p[C.Km5]
    p[C.KimERK] = 1.20e-02
    p[C.KexERK] = 1.80e-02
    p[C.KimpERK] = 1.20e-02
    p[C.KexpERK] = 1.80e-02
    p[C.KimppERK] = 1.10e-02
    p[C.KexppERK] = 1.30e-02
    p[C.V10] = 29.24109258
    p[C.Km10] = 169.0473748
    p[C.n10] = 3.970849295
    p[C.p11] = 0.000126129
    p[C.p12] = 0.007875765
    p[C.p13] = 0.001245747
    p[C.V14] = 5.636949216
    p[C.Km14] = 34180.48
    p[C.V15] = 2.992346912
    p[C.Km15] = 0.001172165
    p[C.p16] = 2.57e-04
    p[C.p17] = 9.63e-05
    p[C.KimDUSP] = 0.024269764
    p[C.KexDUSP] = 0.070467899
    p[C.KimpDUSP] = p[C.KimDUSP]
    p[C.KexpDUSP] = p[C.KexDUSP]
    p[C.V20] = 0.157678678
    p[C.Km20] = 735598.6967
    p[C.V21] = 0.005648117
    p[C.Km21] = 387.8377182
    p[C.p22] = 2.57e-04
    p[C.p23] = 9.63e-05
    p[C.V24] = 0.550346114
    p[C.Km24] = 29516.06587
    p[C.V25] = 10.09063736
    p[C.Km25] = 0.913939859
    p[C.KimRSK] = 0.025925065
    p[C.KexRSK] = 0.129803956
    p[C.V27] = 19.23118154
    p[C.Km27] = 441.5834425
    p[C.V28] = 6.574759504
    p[C.Km28] = 14.99180922
    p[C.V29] = 0.518529841
    p[C.Km29] = 21312.69109
    p[C.V30] = 13.79479021
    p[C.Km30] = 15.04396629
    p[C.V31] = 0.655214248
    p[C.Km31] = 185.9760682
    p[C.n31] = 1.988003164
    p[C.KF31] = 0.013844393
    p[C.nF31] = 2.800340453
    p[C.p32] = 0.003284434
    p[C.p33] = 0.000601234
    p[C.p34] = 7.65E-05
    p[C.V35] = 8.907637012
    p[C.Km35] = 8562.744184
    p[C.V36] = 0.000597315
    p[C.Km36] = 528.552427
    p[C.V37] = 1.745848179
    p[C.Km37] = 0.070379236
    p[C.p38] = 2.57e-04
    p[C.p39] = 9.63e-05
    p[C.KimFOS] = 0.54528521
    p[C.KexFOS] = 0.133249762
    p[C.KimpcFOS] = p[C.KimFOS]
    p[C.KexpcFOS] = p[C.KexFOS]
    p[C.V42] = 0.909968714
    p[C.Km42] = 3992.061328
    p[C.V43] = 0.076717457
    p[C.Km43] = 1157.116021
    p[C.V44] = 0.078344305
    p[C.Km44] = 0.051168202
    p[C.p45] = 2.57e-04
    p[C.p46] = 9.63e-05
    p[C.p47] = 0.001670815
    p[C.m47] = 15.80783969
    p[C.p48] = 0.686020478
    p[C.p49] = 0.314470502
    p[C.m49] = 2.335459127
    p[C.p50] = 26.59483436
    p[C.p51] = 0.01646825
    p[C.m51] = 9.544308421
    p[C.p52] = p[C.p47]
    p[C.m52] = p[C.m47]
    p[C.p53] = p[C.p48]
    p[C.p54] = p[C.p49]
    p[C.m54] = p[C.m49]
    p[C.p55] = p[C.p50]
    p[C.p56] = p[C.p51]
    p[C.m56] = p[C.m51]
    p[C.V57] = 1.026834758
    p[C.Km57] = 0.637490056
    p[C.n57] = 3.584464176
    p[C.p58] = 0.000270488
    p[C.p59] = 0.001443889
    p[C.p60] = 0.002448164
    p[C.p61] = 3.50E-05
    p[C.KimF] = 0.019898797
    p[C.KexF] = 0.396950616
    p[C.p63] = 4.13E-05

    p[C.a] = 218.6276381
    p[C.Vn] = 0.22
    p[C.Vc] = 0.94

    p[C.EGF] = 0.
    p[C.HRG] = 1.
    p[C.no_ligand] = 2.

    return p
end


function initial_values()::Vector{Float64}
    u0::Vector{Float64} = zeros(V.NUM)

    u0[V.ERKc] = 9.60e+02
    u0[V.RSKc] = 3.53e+02
    u0[V.CREBn] = 1.00e+03
    u0[V.Elk1n] = 1.51e+03

    return u0
end
