function MA(k,S)
    return k*S
end


function MAcom2(k,S1,S2)
    return k*S1*S2
end


function MM(k,Km,S)
    return k*S/(Km+S)
end


function Hill2(A,n)
    return A^n
end


function libPulseDelay(t,X,term)
    # signal parameter
    sinput = 1.0 - (1.0 - 0.01)*X
    sbase = 0.01
    slate = 0.05 - (0.05 - 0.01)*X
    tpulse = 0.25
    traise = 0.5
    tdecay = 0.5
    tdelay = 0.0

    if t <= tdelay || term == 1.0
        return sbase
    elseif tdelay < t <= (traise + tdelay)
        return (t - tdelay)*(sinput - sbase)/traise + sbase
    elseif (traise + tdelay) < t <= (tpulse + traise + tdelay)
        return sinput
    elseif (tpulse + traise + tdelay) < t
        return (sinput - slate)*exp(-(t - tpulse - traise - tdelay)/tdecay) + slate
    else
        return 0.0
    end
end


function diffeq(du,u,h,p,t)

    #IKK activity
    IKKa = u[V.IKKp] + u[V.IKKpC] + u[V.IKKppC] + u[V.IKKpp]
    if IKKa == 0.0
        IKKa = 1.0
    end

    ##### time delay
    NFKBnDelay = h(p,t-p[C.delayrnae])[V.NFKBn]

    ###################----CBM module----#######################
    signal = libPulseDelay(t,p[C.X],p[C.term])

    v = Dict{Int64,Float64}()

    v[1] = MA(p[C.kCp0], u[V.C])                                                 # (*basal C phosphorylation*)
    v[2] = signal * MM(p[C.kCpS], p[C.kmCpS], u[V.C])                            # (*signal-dependent C phosphorylation*)
    v[3] = MM(p[C.kCpu], p[C.kmCpu], u[V.Cp])                                    # (* Cp dephosphorylation*)
    v[4] = MA(p[C.kCpB0], u[V.CB])                                               # (*basal CB phosphorylation*)
    v[5] = signal * MM(p[C.kCpBS], p[C.kmCpBS], u[V.CB])                         # (*signal-dependent CB phosphorylation*)
    v[6] = MM(p[C.kCpBu], p[C.kmCpBu], u[V.CpB])                                 # (* CpB dephosphorylation*)
    v[7] = MA(p[C.kCpM0], u[V.CM])                                               # (*basal CM phosphorylation*)
    v[8] = signal * MM(p[C.kCpMS], p[C.kmCpMS], u[V.CM])                         # (*signal-dependent CM phosphorylation*)
    v[9] = MM(p[C.kCpMu], p[C.kmCpMu], u[V.CpM])                                 # (* CpM dephosphorylation*)
    v[10] = MA(p[C.kCpBM0], u[V.CBM])                                            # (*basal CBM phosphorylation*)
    v[11] = signal * MM(p[C.kCpBMS], p[C.kmCpBMS], u[V.CBM])                     # (*signal-dependent CBM phosphorylation*)
    v[12] = MM(p[C.kCpBMu], p[C.kmCpBMu], u[V.CpBM])                             # (* CpBM dephosphorylation*)
    v[13] = MAcom2(p[C.kBaM], u[V.B], u[V.M])                                    # (*association between B and M*)
    v[14] = MA(p[C.kBdM], u[V.BM])                                               # (*disassociation between B and M*)
    v[15] = MAcom2(p[C.kCaB], u[V.C], u[V.B])                                    # (*association between C and B*)
    v[16] = MA(p[C.kCdB], u[V.CB])                                               # (*disassociation between C and B*)
    v[17] = MAcom2(p[C.kCpaB], u[V.Cp], u[V.B])                                  # (*association between Cp and B*)
    v[18] = MA(p[C.kCpdB], u[V.CpB])                                             # (*disassociation between Cp and B*)
    v[19] = MAcom2(p[C.kCaM], u[V.C], u[V.M])                                    # (*association between C and M*)
    v[20] = MA(p[C.kCdM], u[V.CM])                                               # (*disassociation between C and M*)
    v[21] = MAcom2(p[C.kCpaM], u[V.Cp], u[V.M])                                  # (*association between Cp and M*)
    v[22] = MA(p[C.kCpdM], u[V.CpM])                                             # (*disassociation between Cp and M*)
    v[23] = MAcom2(p[C.kCBaM], u[V.CB], u[V.M])                                  # (*association between CB and M*)
    v[24] = MA(p[C.kCBdM], u[V.CBM])                                             # (*disassociation between CB and M*)
    v[25] = MAcom2(p[C.kCpBaM], u[V.CpB], u[V.M])                                # (*association between CpB and M*)
    v[26] = MA(p[C.kCpBdM], u[V.CpBM])                                           # (*disassociation between CpB and M*)
    v[27] = MAcom2(p[C.kCMaB], u[V.CM], u[V.B])                                  # (*association between CM and B*)
    v[28] = MA(p[C.kCMdB], u[V.CBM])                                             # (*disassociation between CM and B*)
    v[29] = MAcom2(p[C.kCpMaB], u[V.CpM], u[V.B])                                # (*association between CpM and B*)
    v[30] = MA(p[C.kCpMdB], u[V.CpBM])                                           # (*disassociation between CpM and B*)
    v[31] = MAcom2(p[C.kCaBM], u[V.C], u[V.BM])                                  # (*association between C and BM*)
    v[32] = MA(p[C.kCdBM], u[V.CBM])                                             # (*disassociation between C and BM*)
    v[33] = MAcom2(p[C.kCpaBM], u[V.Cp], u[V.BM])                                # (*association between Cp and BM*)
    v[34] = MA(p[C.kCpdBM], u[V.CpBM])                                           # (*disassociation between Cp and BM*)
    ###################----CBM module----#######################

    ###################----TAK1 module----######################
    v[35] = MA(p[C.kTp0], u[V.TAK1])                                             # (*basal TAK1 phosphorylation*)
    v[36] = u[V.Cp] * MM(p[C.kCTpS], p[C.kmCTpS], u[V.TAK1])                     # (*Cp-dependent TAK1 phosphorylation*)
    v[37] = u[V.CpB] * MM(p[C.kCBTpS], p[C.kmCBTpS], u[V.TAK1])                  # (*CpB-dependent TAK1 phosphorylation*)
    v[38] = u[V.CpM] * MM(p[C.kCMTpS], p[C.kmCMTpS], u[V.TAK1])                  # (*CpM-dependent TAK1 phosphorylation*)
    v[39] = u[V.CpBM] * MM(p[C.kCBMTpS], p[C.kmCBMTpS], u[V.TAK1])               # (*CpBM-dependent TAK1 phosphorylation*)
    v[40] = u[V.IKKp] * MM(p[C.kTpIKK1], p[C.kmTpIKK1], u[V.TAK1])               # (*IKKp-dependent TAK1 phosphorylation*)
    v[41] = MM(p[C.kTpu], p[C.kmTpu], u[V.TAK1p])                                # (*TAK1 dephosphorylation*)
    v[42] = MA(p[C.kTpC0], u[V.TAK1C])                                           # (*basal TAK1C phosphorylation*)
    v[43] = u[V.Cp] * MM(p[C.kCTpCS], p[C.kmCTpCS], u[V.TAK1C])                  # (*Cp-dependent TAK1C phosphorylation*)
    v[44] = u[V.CpB] * MM(p[C.kCBTpCS], p[C.kmCBTpCS], u[V.TAK1C])               # (*CpB-dependent TAK1C phosphorylation*)
    v[45] = u[V.CpM] * MM(p[C.kCMTpCS], p[C.kmCMTpCS], u[V.TAK1C])               # (*CpM-dependent TAK1C phosphorylation*)
    v[46] = u[V.CpBM] * MM(p[C.kCBMTpCS], p[C.kmCBMTpCS], u[V.TAK1C])            # (*CpBM-dependent TAK1C phosphorylation*)
    v[47] = u[V.IKKp] * MM(p[C.kTpCIKK1], p[C.kmTpCIKK1], u[V.TAK1C])            # (*IKKp-dependent TAK1C phosphorylation*)
    v[48] = u[V.IKKpC] * MM(p[C.kTpCIKK2], p[C.kmTpCIKK2], u[V.TAK1C])           # (*IKKpC-dependent TAK1C phosphorylation*)
    v[49] = u[V.IKKppC] * MM(p[C.kTpCIKK3], p[C.kmTpCIKK3], u[V.TAK1C])          # (*IKKppC-dependent TAK1C phosphorylation*)
    v[50] = MM(p[C.kTpCu], p[C.kmTpCu], u[V.TAK1pC])                             # (*basal TAK1C dephosphorylation*)
    v[51] = MAcom2(p[C.kTaB], (u[V.CB] + u[V.CpB]), u[V.TAK1])                   # (*association between TAK1 and CB*)
    v[52] = MAcom2(p[C.kTaM], (u[V.CM] + u[V.CpM]), u[V.TAK1])                   # (*association between TAK1 and CM*)
    v[53] = MAcom2(p[C.kTaBM], (u[V.CBM] + u[V.CpBM]), u[V.TAK1])                # (*association between TAK1 and CBM*)
    v[54] = MA(p[C.kTCd], u[V.TAK1C])                                            # (*dissociation between TAK1 and CBM*)
    v[55] = MAcom2(p[C.kTpaB], (u[V.CB] + u[V.CpB]), u[V.TAK1p])                 # (*association between TAK1p and CBM*)
    v[56] = MAcom2(p[C.kTpaM], (u[V.CM] + u[V.CpM]), u[V.TAK1p])                 # (*association between TAK1p and CBM*)
    v[57] = MAcom2(p[C.kTpaBM], (u[V.CBM] + u[V.CpBM]), u[V.TAK1p])              # (*association between TAK1p and CBM*)
    v[58] = MA(p[C.kTpCd], u[V.TAK1pC])                                          # (*dissociation between TAK1 and CBM*)
    ###################----TAK1 module----######################

    ###################----IKK module----#######################
    v[59] = MA(p[C.kIp0], u[V.IKK])                                              # (*basal IKK phosphorylation*)
    v[60] = u[V.TAK1p] * MM(p[C.kIpTAKp], p[C.kmIpTAKp], u[V.IKK])               # (*TAK1p-dependent IKK phosphorylation*)
    v[61] = u[V.TAK1pC] * MM(p[C.kIpTAKpC], p[C.kmIpTAKpC], u[V.IKK])            # (*TAKpC-dependent IKK phosphorylation*)
    v[62] = MA(p[C.kICp0], u[V.IKKC])                                            # (*basal IKKC phosphorylation*)
    v[63] = u[V.TAK1p] * MM(p[C.kICpTAKp], p[C.kmICpTAKp], u[V.IKKC])            # (*TAK1p-dependent IKKC phosphorylation*)
    v[64] = u[V.TAK1pC] * MM(p[C.kICpTAKpC], p[C.kmICpTAKpC], u[V.IKKC])         # (*TAK1pC-dependent IKKC phosphorylation*)
    v[65] = MM(p[C.kIpu], p[C.kmIpu], u[V.IKKp])                                 # (*IKKp dephosphorylation*)
    v[66] = MM(p[C.kIpCu], p[C.kmIpCu], u[V.IKKpC])                              # (*IKKpC dephosphorylation*)
    v[67] = MAcom2(p[C.kIaB], (u[V.CB] + u[V.CpB]), u[V.IKK])                    # (*association between IKK and CB*)
    v[68] = MAcom2(p[C.kIaM], (u[V.CM] + u[V.CpM]), u[V.IKK])                    # (*association between IKK and CM*)
    v[69] = MAcom2(p[C.kIaBM], (u[V.CBM] + u[V.CpBM]), u[V.IKK])                 # (*association between IKK and CBM*)
    v[70] = MA(p[C.kICd], u[V.IKKC])                                             # (*dissociation between IKK and CBM*)
    v[71] = MAcom2(p[C.kIpaB], (u[V.CB] + u[V.CpB]), u[V.IKKp])                  # (*association between IKKp and CB*)
    v[72] = MAcom2(p[C.kIpaM], (u[V.CM] + u[V.CpM]), u[V.IKKp])                  # (*association between IKKp and CM*)
    v[73] = MAcom2(p[C.kIpaBM], (u[V.CBM] + u[V.CpBM]), u[V.IKKp])               # (*association between IKKp and CBM*)
    v[74] = MA(p[C.kIpCd], u[V.IKKpC])                                           # (*dissociation between IKKp and CBM*)
    v[75] = MM(p[C.kIpCfaIKKpC], p[C.kmIpCfaIKKpC], u[V.IKKpC])                  # (*basal IKKpC phosphorylation*)
    v[76] = u[V.IKKppC] * MM(p[C.kIpCfaIKKppC], p[C.kmIpCfaIKKppC], u[V.IKKpC])  # (*IKKppC-dependent IKKpC phosphorylation*)
    v[77] = MM(p[C.kIppCu], p[C.kmIppCu], u[V.IKKppC])                           # (*basal IKKppC dephosphorylation*)
    v[78] = MA(p[C.kIppCd], u[V.IKKppC])                                         # (*dissociation between IKKpp and CBM (IKKppC->IKKpp*)
    v[79] = MM(p[C.kIpphf], p[C.kmIpphf], u[V.IKKpp])                            # (*basal IKKpp inactivation (IKKpp->IKKi)*)
    v[80] = MM(p[C.kIppChf], p[C.kmIppChf], u[V.IKKppC])                         # (*basal IKKppC inactivation (IKKppC->IKKi)*)
    v[81] = MM(p[C.kIir], p[C.kmIir], u[V.IKKi])                                 # (* recycling (IKKi->IKK)*)
    v[82] = u[V.A20c] * MM(p[C.kIpA20], p[C.kmIpA20], u[V.IKKp])                 # (* A20-dependent IKKp phosphorylation inhibition *)
    v[83] = u[V.A20c] * MM(p[C.kIpCA20], p[C.kmIpCA20], u[V.IKKpC])              # (* A20-dependent IKKpC phosphorylation inhibition *)
    ###################----IKK module----#######################

    ########----IKKNFkBIkB module----################
    v[84] = MAcom2(p[C.kassanfkbikk], IKKa, u[V.NFKBIKBac])                      # (*association between IKKa and NFkBIkBac*)
    v[85] = MAcom2(p[C.kassbnfkbikk], IKKa, u[V.NFKBIKBbc])                      # (*association between IKKa and NFkBIkBbc*)
    v[86] = MAcom2(p[C.kassenfkbikk], IKKa, u[V.NFKBIKBec])                      # (*association between IKKa and NFkBIkBec*)
    v[87] = MA(p[C.kdisanfkbikk], u[V.IKKNFKBIKBac])                             # (*dissociation between IKKa and NFkBIkBac*)
    v[88] = MA(p[C.kdisbnfkbikk], u[V.IKKNFKBIKBbc])                             # (*dissociation between IKKa and NFkBIkBbc*)
    v[89] = MA(p[C.kdisenfkbikk], u[V.IKKNFKBIKBec])                             # (*dissociation between IKKa and NFkBIkBec*)
    v[90] = MA(p[C.kdegboundaIKK], u[V.IKKNFKBIKBac])                            # (*dissociation between IKKa and NFkBc and degradation of IkBac*)
    v[91] = MA(p[C.kdegboundbIKK], u[V.IKKNFKBIKBbc])                            # (*dissociation between IKKa and NFkBc and degradation of IkBbc*)
    v[92] = MA(p[C.kdegboundeIKK], u[V.IKKNFKBIKBec])                            # (*dissociation between IKKa and NFkBc and degradation of IkBec*)
    v[93] = MAcom2(p[C.kassaikk], u[V.IKBac], IKKa)                              # (*association between IKKa and IkBac*)
    v[94] = MAcom2(p[C.kassbikk], u[V.IKBbc], IKKa)                              # (*association between IKKa and IkBbc*)
    v[95] = MAcom2(p[C.kasseikk], u[V.IKBec], IKKa)                              # (*association between IKKa and IkBec*)
    v[96] = MA(p[C.kdisaikk], u[V.IKKIKBac])                                     # (*dissociation between IKKa and IkBac*)
    v[97] = MA(p[C.kdisbikk], u[V.IKKIKBbc])                                     # (*dissociation between IKKa and IkBbc*)
    v[98] = MA(p[C.kdiseikk], u[V.IKKIKBec])                                     # (*dissociation between IKKa and IkBec*)
    v[99] = MA(p[C.kdegfreeaIKK], u[V.IKKIKBac])                                 # (*dissociation of IKKa and degradation of IkBac*)
    v[100] = MA(p[C.kdegfreebIKK], u[V.IKKIKBbc])                                # (*dissociation of IKKa and degradation of IkBbc*)
    v[101] = MA(p[C.kdegfreeeIKK], u[V.IKKIKBec])                                # (*dissociation of IKKa and degradation of IkBec*)
    ########----IKKNFkBIkB module----################

    ###################----NFkB module----######################
    v[102] = MAcom2(p[C.kassa], u[V.IKBac], u[V.NFKBc])                          # (*association between NFkBc and IkBac*)
    v[103] = MA(p[C.kdisa], u[V.NFKBIKBac])                                      # (*dissociation between NFkBc and IkBac*)
    v[104] = MA(p[C.kdegbounda], u[V.NFKBIKBac])                                 # (*dissociation of NFkBc and degradation of IkBac*)
    v[105] = MAcom2(p[C.kassaikknfkb], u[V.IKKIKBac], u[V.NFKBc])                # (*association between NFkBc and IKKaIkBac*)
    v[106] = MA(p[C.kdisaikknfkb], u[V.IKKNFKBIKBac])                            # (*dissociation between NFkBc and IKKaIkBac*)
    v[107] = MAcom2(p[C.kassa], u[V.IKBan], u[V.NFKBn])                          # (*association between NFkBn and IkBan*)
    v[108] = MA(p[C.kdisa], u[V.NFKBIKBan])                                      # (*dissociation between NFkBn and IkBan*)
    v[109] = MA(p[C.kdegbounda], u[V.NFKBIKBan])                                 # (*dissociation of NFkBn and degradation of IkBan*)

    v[110] = MAcom2(p[C.kassb], u[V.IKBbc], u[V.NFKBc])                          # (*association between NFkBc and IkBbc*)
    v[111] = MA(p[C.kdisb], u[V.NFKBIKBbc])                                      # (*dissociation between NFkBc and IkBbc*)
    v[112] = MA(p[C.kdegboundb], u[V.NFKBIKBbc])                                 # (*dissociation of NFkBc and degradation of IkBbc*)
    v[113] = MAcom2(p[C.kassbikknfkb], u[V.IKKIKBbc], u[V.NFKBc])                # (*association between NFkBc and IKKaIkBbc*)
    v[114] = MA(p[C.kdisbikknfkb], u[V.IKKNFKBIKBbc])                            # (*dissociation between NFkBc and IKKaIkBbc*)
    v[115] = MAcom2(p[C.kassb], u[V.IKBbn], u[V.NFKBn])                          # (*association between NFkBn and IkBbn*)
    v[116] = MA(p[C.kdisb], u[V.NFKBIKBbn])                                      # (*dissociation between NFkBn and IkBbn*)
    v[117] = MA(p[C.kdegboundb], u[V.NFKBIKBbn])                                 # (*dissociation of NFkBn and degradation of IkBbn*)

    v[118] = MAcom2(p[C.kasse], u[V.IKBec], u[V.NFKBc])                          # (*association between NFkBc and IkBec*)
    v[119] = MA(p[C.kdise], u[V.NFKBIKBec])                                      # (*dissociation between NFkBc and IkBec*)
    v[120] = MA(p[C.kdegbounde], u[V.NFKBIKBec])                                 # (*dissociation of NFkBc and degradation of IkBec*)
    v[121] = MAcom2(p[C.kasseikknfkb], u[V.IKKIKBec], u[V.NFKBc])                # (*association between NFkBc and IKKaIkBec*)
    v[122] = MA(p[C.kdiseikknfkb], u[V.IKKNFKBIKBec])                            # (*dissociation between NFkBc and IKKaIkBec*)
    v[123] = MAcom2(p[C.kasse], u[V.IKBen], u[V.NFKBn])                          # (*association between NFkBn and IkBen*)
    v[124] = MA(p[C.kdise], u[V.NFKBIKBen])                                      # (*dissociation between NFkBn and IkBen*)
    v[125] = MA(p[C.kdegbounde], u[V.NFKBIKBen])                                 # (*dissociation of NFkBn and degradation of IkBen*)

    v[126] = MA(p[C.kshutboundikbain], u[V.NFKBIKBac])                           # (*transport NFkBIkBac into nucleos*)
    v[127] = MA(p[C.kshutboundikbaout], u[V.NFKBIKBan])                          # (*transport NFkBIkBan into cytoplasm*)
    v[128] = MA(p[C.kshutboundikbbin], u[V.NFKBIKBbc])                           # (*transport NFkBIkBbc into nucleos*)
    v[129] = MA(p[C.kshutboundikbbout], u[V.NFKBIKBbn])                          # (*transport NFkBIkBbn into cytoplasm*)
    v[130] = MA(p[C.kshutboundikbein], u[V.NFKBIKBec])                           # (*transport NFkBIkBec into nucleos*)
    v[131] = MA(p[C.kshutboundikbeout], u[V.NFKBIKBen])                          # (*transport NFkBIkBen into cytoplasm*)
    v[132] = MA(p[C.kshutfreeikbain], u[V.IKBac])                                # (*transport IkBac into nucleos*)
    v[133] = MA(p[C.kshutfreeikbaout], u[V.IKBan])                               # (*transport IkBan into cytoplasm*)
    v[134] = MA(p[C.kshutfreeikbbin], u[V.IKBbc])                                # (*transport IkBbc into nucleos*)
    v[135] = MA(p[C.kshutfreeikbbout], u[V.IKBbn])                               # (*transport IkBbn into cytoplasm*)
    v[136] = MA(p[C.kshutfreeikbein], u[V.IKBec])                                # (*transport IkBec into nucleos*)
    v[137] = MA(p[C.kshutfreeikbeout], u[V.IKBen])                               # (*transport IkBen into cytoplasm*)
    v[138] = MA(p[C.kshutnfkbin], u[V.NFKBc])                                    # (*transport NFkBc into nucleos*)
    v[139] = MA(p[C.kshutnfkbout], u[V.NFKBn])                                   # (*transport NFkBn into cytoplasm*)

    v[140] = p[C.k0mrnaikba]                                                     # (*basal mRNA(IkBa) transcription*)
    v[141] = p[C.kprodmrnaikba] * Hill2(u[V.NFKBn], p[C.khillprodmrnaikba])      # (*NFkB-induced IkBa transcription*)
    v[142] = MA(p[C.kdegmrnaikba], u[V.mRNAac])                                  # (*mRNA(IkBa) degradation*)
    v[143] = MA(p[C.kpikba], u[V.mRNAac])                                        # (*mRNA(IkBa) translation*)
    v[144] = MA(p[C.kdegfreea], u[V.IKBac])                                      # (*IkBac degradation*)
    v[145] = MA(p[C.kdegfreea], u[V.IKBan])                                      # (*IkBan degradation*)

    v[146] = p[C.k0mrnaikbb]                                                     # (*basal mRNA(IkBb) transcription*)
    v[147] = MA(p[C.kdegmrnaikbb], u[V.mRNAbc])                                  # (*mRNA(IkBb) degradation*)
    v[148] = MA(p[C.kpikbb], u[V.mRNAbc])                                        # (*mRNA(IkBb) translation*)
    v[149] = MA(p[C.kdegfreeb], u[V.IKBbc])                                      # (*IkBbc degradation*)
    v[150] = MA(p[C.kdegfreeb], u[V.IKBbn])                                      # (*IkBbn degradation*)

    v[151] = p[C.k0mrnaikbe]                                                     # (*basal mRNA(IkBe) transcription*)
    v[152] = p[C.kprodmrnaikbe] * Hill2(NFKBnDelay, p[C.khillprodmrnaikbe])      # (*NFkB-induced IkBe transcription*)
    v[153] = MA(p[C.kdegmrnaikbe], u[V.mRNAec])                                  # (*mRNA(IkBe) degradation*)
    v[154] = MA(p[C.kpikbe], u[V.mRNAec])                                        # (*mRNA(IkBe) translation*)
    v[155] = MA(p[C.kdegfreee], u[V.IKBec])                                      # (*IkBec degradation*)
    v[156] = MA(p[C.kdegfreee], u[V.IKBen])                                      # (*IkBen degradation*)

    v[157] = p[C.k0mrnaa20]                                                      # (*basal mRNA(A20) transcription*)
    v[158] = p[C.kprodmrnaa20] * Hill2(u[V.NFKBn], p[C.khillprodmrnaa20])        # (*NFkB-induced A20 transcription*)
    v[159] = MA(p[C.kdegmrnaa20], u[V.mRNAa20c])                                 # (*mRNA(A20) degradation*)
    v[160] = MA(p[C.kpa20], u[V.mRNAa20c])                                       # (*mRNA(A20) translation*)
    v[161] = MA(p[C.kdega20], u[V.A20c])                                         # (*A20 degradation*)

    IKKaNFKBIKB = - v[93] - v[94] - v[95] + v[99] + v[96] + v[100] + v[97] + v[101] + v[98] + v[90] + v[87] +
                    v[91] + v[88] + v[92] + v[89] - v[84] - v[85] - v[86]

    ###################----CBM module----######################
    du[V.B] = (- v[15] + v[16]) + (- v[17] + v[18]) + (- v[27] + v[28]) + (- v[29] + v[30]) + (- v[13] + v[14])
    du[V.M] = (- v[19] + v[20]) + (- v[21] + v[22]) + (- v[23] + v[24]) + (- v[25] + v[26]) + (- v[13] + v[14])
    du[V.BM] = ( v[13] - v[14]) + (- v[31] + v[32]) + (- v[33] + v[34])
    du[V.C] = (- v[15] + v[16]) + (- v[19] + v[20]) + (- v[31] + v[32]) + ( v[3] -(v[1] + v[2]) )
    du[V.CB] = (v[15] - v[16]) + (- v[23] + v[24]) + ( v[6] -(v[4] + v[5]) )
    du[V.CM] = (v[19] - v[20]) + (- v[27] + v[28]) + ( v[9] -(v[7] + v[8]) )
    du[V.Cp] = (- v[17] + v[18]) + (- v[21] + v[22]) + (- v[33] + v[34]) + (- v[3] +(v[1] + v[2]) )
    du[V.CpB] = (v[17] - v[18]) + (- v[25] + v[26]) + (- v[6] +(v[4] + v[5]) )
    du[V.CpM] = (v[21] - v[22]) + (- v[29] + v[30]) + (- v[9] +(v[7] + v[8]) )
    du[V.CBM] = (v[27] - v[28]) + (v[23] - v[24]) + ( v[31] - v[32]) + ( v[12] -(v[10] + v[11]) )
    du[V.CpBM] = (v[29] - v[30]) + (v[25] - v[26]) + ( v[33] - v[34]) + ( - v[12] +(v[10] + v[11]) )
    ###################----CBM module----######################

    ###################----TAK1 module---######################
    du[V.TAK1] =  -(v[35] + v[36] + v[37] + v[38] + v[39] + v[40])  + v[41] - (v[51] + v[52] + v[53]) + v[54]
    du[V.TAK1p] =  (v[35] + v[36] + v[37] + v[38] + v[39] + v[40]) - v[41] - (v[55] + v[56] + v[57]) + v[58]
    du[V.TAK1C] = -(v[42] + v[43] + v[44] + v[45] + v[46] + v[47] + v[48] + v[49] ) + v[50] + (v[51] + v[52] + v[53]) - v[54]
    du[V.TAK1pC] =  (v[42] + v[43] + v[44] + v[45] + v[46] + v[47] + v[48] + v[49] ) - v[50] + (v[55] + v[56] + v[57]) - v[58]
    ###################----TAK1 module---######################

    ###################----IKK module----######################
    du[V.IKK] =  -(v[59] + v[60] + v[61]) + v[65] + v[81] - (v[67] + v[68] + v[69]) + v[70] + v[82]
    du[V.IKKC] = -(v[62] + v[63] + v[64]) + v[66] + (v[67] + v[68] + v[69]) - v[70] + v[83]
    du[V.IKKp] = (v[59] + v[60] + v[61]) - v[65] - (v[71] + v[72] + v[73]) + v[74] - v[82] + u[V.IKKp]*IKKaNFKBIKB/IKKa
    du[V.IKKpC] = (v[71] + v[72] + v[73]) - v[74] - ( v[75] + v[76]) + v[77]+ (v[62] + v[63] + v[64]) - v[66] - v[83] + u[V.IKKpC]*IKKaNFKBIKB/IKKa
    du[V.IKKppC] =  ( v[75] + v[76]) - v[77] - v[78] - v[80] + u[V.IKKppC]*IKKaNFKBIKB/IKKa
    du[V.IKKpp] = v[78] - v[79] + u[V.IKKpp]*IKKaNFKBIKB/IKKa
    du[V.IKKi] = v[79] + v[80] - v[81]
    ###################----IKK module----######################

    ###################----NFkB module----######################
    du[V.NFKBc] = v[90] + v[91] + v[92] - v[102] + v[103]+ v[104] - v[105] + v[106] - v[110] + v[111] +
                v[112] - v[113] + v[114] - v[118] + v[119] + v[120] - v[121] + v[122] - v[138] + v[139]
    du[V.NFKBn] = - v[107] + v[108] + v[109] - v[115] + v[116] + v[117] - v[123] + v[124] + v[125] + v[138] - v[139]
    du[V.IKBac] = - v[93] + v[96] - v[102] + v[103] - v[132] + v[133] + v[143] - v[144]
    du[V.IKBbc] = - v[94] + v[97] - v[110] + v[111] - v[134] + v[135] + v[148] - v[149]
    du[V.IKBec] = - v[95] + v[98] - v[118] + v[119] - v[136] + v[137] + v[154] - v[155]
    du[V.IKBan] = - v[107] + v[108] + v[132] - v[133] - v[145]
    du[V.IKBbn] = - v[115] + v[116] + v[134] - v[135] - v[150]
    du[V.IKBen] = - v[123] + v[124] + v[136] - v[137] - v[156]

    du[V.NFKBIKBac] = - v[84] + v[87] + v[102] - v[103] - v[104] - v[126] + v[127]
    du[V.NFKBIKBbc] = - v[85] + v[88] + v[110] - v[111] - v[112] - v[128] + v[129]
    du[V.NFKBIKBec] = - v[86] + v[89] + v[118] - v[119] - v[120] - v[130] + v[131]
    du[V.NFKBIKBan] = v[107] - v[108] - v[109] + v[126] - v[127]
    du[V.NFKBIKBbn] = v[115] - v[116] - v[117] + v[128] - v[129]
    du[V.NFKBIKBen] = v[123] - v[124] - v[125] + v[130] - v[131]
    du[V.IKKIKBac] = v[93] - v[96] - v[99] - v[105] + v[106]
    du[V.IKKIKBbc] = v[94] - v[97] - v[100] - v[113] + v[114]
    du[V.IKKIKBec] = v[95] - v[98] - v[101] - v[121] + v[122]
    du[V.IKKNFKBIKBac] =  v[84] - v[87] - v[90] + v[105] - v[106]
    du[V.IKKNFKBIKBbc] =  v[85] - v[88] - v[91] + v[113] - v[114]
    du[V.IKKNFKBIKBec] =  v[86] - v[89] - v[92] + v[121] - v[122]

    du[V.mRNAac] = v[140] + v[141] - v[142]
    du[V.mRNAbc] = v[146] - v[147]
    du[V.mRNAec] = v[151] + v[152] - v[153]
    du[V.mRNAa20c] = v[157] + v[158] - v[159]
    du[V.A20c] = v[160] - v[161]
    ###################----NFkB module----######################
end


function param_values()

    p = zeros(C.NUM)

    p[C.kCp0] = 1.000000e-04
    p[C.kCpB0] = 1.000000e-04
    p[C.kCpM0] = 1.000000e-04
    p[C.kCpBM0] = 1.000000e-04
    p[C.kCpS] = 7.000000e+00
    p[C.kCpBS] = 2.000000e+01
    p[C.kCpMS] = 1.700000e+01
    p[C.kCpBMS] = 1.000000e+00
    p[C.kmCpS] = 1.000000e+00
    p[C.kmCpBS] = 1.000000e+00
    p[C.kmCpMS] = 1.000000e+00
    p[C.kmCpBMS] = 5.000000e+00
    p[C.kCpu] = 7.500000e-01
    p[C.kCpBu] = 1.950000e+00
    p[C.kCpMu] = 3.750000e+00
    p[C.kCpBMu] = 2.500000e-01
    p[C.kmCpu] = 1.000000e+00
    p[C.kmCpBu] = 1.000000e+00
    p[C.kmCpMu] = 1.000000e+00
    p[C.kmCpBMu] = 2.500000e+00

    p[C.kBaM] = 5.000000e-01
    p[C.kCaB] = 5.000000e-01
    p[C.kCaM] = 5.000000e-01
    p[C.kCBaM] = 1.500000e+00
    p[C.kCMaB] = 5.000000e-01
    p[C.kCaBM] = 5.000000e-01
    p[C.kCpaBM] = 1.885000e+00
    p[C.kCpaB] = 7.885000e+00
    p[C.kCpaM] = 4.585000e+00
    p[C.kCpBaM] = 2.500000e+00
    p[C.kCpMaB] = 1.500000e+00
    p[C.kBdM] = 2.500000e+00
    p[C.kCdB] = 1.000000e+00
    p[C.kCdM] = 2.500000e+00
    p[C.kCBdM] = 1.500000e+00
    p[C.kCMdB] = 2.500000e+00
    p[C.kCpdB] = 3.000000e-01
    p[C.kCpdM] = 3.000000e-01
    p[C.kCpBdM] = 1.500000e-01
    p[C.kCpMdB] = 4.500000e-01
    p[C.kCdBM] = 2.000000e-01
    p[C.kCpdBM] = 2.500000e-01

    p[C.kTp0] = 2.700000e-01
    p[C.kCTpS] = 3.200000e-01
    p[C.kCBTpS] = 1.320000e+00
    p[C.kCMTpS] = 1.320000e+00
    p[C.kCBMTpS] = 1.320000e+00
    p[C.kmCTpS] = 2.100000e+00
    p[C.kmCBTpS] = 2.100000e+00
    p[C.kmCMTpS] = 2.100000e+00
    p[C.kmCBMTpS] = 2.100000e+00
    p[C.kTpIKK1] = 1.670000e+01
    p[C.kmTpIKK1] = 7.800000e-01
    p[C.kTpu] = 3.400000e+00
    p[C.kmTpu] = 1.790000e+00
    p[C.kTpC0] = 3.300000e-01
    p[C.kCTpCS] = 5.500000e+00
    p[C.kCBTpCS] = 1.750000e+01
    p[C.kCMTpCS] = 1.650000e+01
    p[C.kCBMTpCS] = 2.450000e+01
    p[C.kmCTpCS] = 5.200000e-01
    p[C.kmCBTpCS] = 5.200000e-01
    p[C.kmCMTpCS] = 5.200000e-01
    p[C.kmCBMTpCS] = 5.200000e-01
    p[C.kTpCIKK1] = 1.000000e+01
    p[C.kmTpCIKK1] = 7.800000e-01
    p[C.kTpCIKK2] = 1.710000e+01
    p[C.kmTpCIKK2] = 1.300000e-01
    p[C.kTpCIKK3] = 8.000000e+01
    p[C.kmTpCIKK3] = 7.000000e-01
    p[C.kTpCu] = 1.533000e+01
    p[C.kmTpCu] = 8.000000e-01
    p[C.kTaB] = 1.790000e+02
    p[C.kTaM] = 1.690000e+02
    p[C.kTaBM] = 3.690000e+02
    p[C.kTpaB] = 3.790000e+00
    p[C.kTpaM] = 3.790000e+00
    p[C.kTpaBM] = 3.790000e+00
    p[C.kTCd] = 4.000000e+02
    p[C.kTpCd] = 3.210000e+00

    p[C.kIp0] = 5.500000e-02
    p[C.kIpTAKp] = 1.930000e+00
    p[C.kmIpTAKp] = 4.000000e-01
    p[C.kIpTAKpC] = 3.800000e+00
    p[C.kmIpTAKpC] = 5.800000e-01
    p[C.kIpu] = 1.500000e+01
    p[C.kmIpu] = 6.700000e-01
    p[C.kICp0] = 5.100000e-02
    p[C.kICpTAKp] = 8.640000e+00
    p[C.kmICpTAKp] = 7.100000e+00
    p[C.kICpTAKpC] = 1.230000e+01
    p[C.kmICpTAKpC] = 7.177000e+00
    p[C.kIpCu] = 1.250000e+01
    p[C.kmIpCu] = 5.410000e+00
    p[C.kIaB] = 3.030000e+01
    p[C.kIaM] = 2.230000e+01
    p[C.kIaBM] = 5.030000e+01
    p[C.kIpaB] = 3.130000e+01
    p[C.kIpaM] = 2.330000e+01
    p[C.kIpaBM] = 6.030000e+01
    p[C.kICd] = 3.280000e+01
    p[C.kIpCd] = 1.320000e+00
    p[C.kIpCfaIKKpC] = 1.000000e-05
    p[C.kmIpCfaIKKpC] = 2.560000e+00
    p[C.kIpCfaIKKppC] = 1.085700e+02
    p[C.kmIpCfaIKKppC] = 2.010000e+00
    p[C.kIppCu] = 1.530000e+00
    p[C.kmIppCu] = 2.000000e-01
    p[C.kIppCd] = 4.500000e-01
    p[C.kIppChf] = 1.000000e-01
    p[C.kmIppChf] = 5.000000e-01
    p[C.kIpphf] = 9.000000e-03
    p[C.kmIpphf] = 1.520000e+00
    p[C.kIir] = 1.440000e+00
    p[C.kmIir] = 3.440000e+00
    p[C.kIpA20] = 2.000000e+01
    p[C.kmIpA20] = 1.000000e-03
    p[C.kIpCA20] = 2.000000e+01
    p[C.kmIpCA20] = 1.000000e-03


    p[C.kassaikk] = 1.631051e-03
    p[C.kassaikknfkb] = 1.305073e+01
    p[C.kassanfkbikk] = 4.812297e-02
    p[C.kdisaikk] = 1.375438e-01
    p[C.kdisaikknfkb] = 2.296074e-05
    p[C.kdisanfkbikk] = 3.030804e-02
    p[C.kassa] = 3.655315e+01
    p[C.kdisa] = 3.495116e-04
    p[C.kdegbounda] = 2.089360e-03
    p[C.kdegboundaIKK] = 1.774837e-01
    p[C.kdegfreea] = 2.010515e-02
    p[C.kdegfreeaIKK] = 2.265038e-03
    p[C.kshutboundikbain] = 4.102701e-01
    p[C.kshutboundikbaout] = 1.632918e-01
    p[C.kshutfreeikbain] = 1.011242e-01
    p[C.kshutfreeikbaout] = 1.158262e-03
    p[C.k0mrnaikba] = 3.450000e-05
    p[C.kdegmrnaikba] = 8.326217e-02
    p[C.khillprodmrnaikba] = 1.620503e+00
    p[C.kpikba] = 8.484488e+00
    p[C.kprodmrnaikba] = 1.290500e-01

    p[C.kassbikk] = 5.609029e-03
    p[C.kassbikknfkb] = 1.423728e+02
    p[C.kassbnfkbikk] = 1.541680e-02
    p[C.kdisbikk] = 2.771705e-02
    p[C.kdisbikknfkb] = 1.742842e-03
    p[C.kdisbnfkbikk] = 1.068148e-01
    p[C.kassb] = 2.437209e+01
    p[C.kdegboundb] = 2.804468e-04
    p[C.kdegboundbIKK] = 3.645044e-01
    p[C.kdegfreeb] = 1.486087e-01
    p[C.kdegfreebIKK] = 2.818781e-04
    p[C.kdisb] = 1.856372e-04
    p[C.kshutboundikbbin] = 2.155431e-01
    p[C.kshutboundikbbout] = 7.783843e-01
    p[C.kshutfreeikbbin] = 1.505364e-03
    p[C.kshutfreeikbbout] = 8.350086e-04
    p[C.k0mrnaikbb] = 9.336099e-05
    p[C.kdegmrnaikbb] = 4.101465e-03
    p[C.khillprodmrnaikbb] = 0.000000e+00
    p[C.kpikbb] = 4.171048e-02
    p[C.kprodmrnaikbb] = 0.000000e+00
    p[C.kasseikk] = 4.423523e-03
    p[C.kasseikknfkb] = 1.485352e+02
    p[C.kassenfkbikk] = 1.531672e-02
    p[C.kdiseikk] = 3.827432e-01
    p[C.kdiseikknfkb] = 2.124084e-04
    p[C.kdisenfkbikk] = 2.691468e-01

    p[C.kasse] = 2.413478e+02
    p[C.kdegbounde] = 1.712910e-02
    p[C.kdegboundeIKK] = 6.934859e-02
    p[C.kdegfreee] = 2.201754e-01
    p[C.kdegfreeeIKK] = 4.819481e-05
    p[C.kdise] = 2.215707e-02
    p[C.kshutboundikbein] = 8.666062e-03
    p[C.kshutboundikbeout] = 4.006261e-02
    p[C.kshutfreeikbein] = 1.364442e-03
    p[C.kshutfreeikbeout] = 1.014252e-04
    p[C.k0mrnaikbe] = 1.232627e-04
    p[C.kdegmrnaikbe] = 9.773320e-03
    p[C.khillprodmrnaikbe] = 1.944545e+00
    p[C.kpikbe] = 2.384772e-01
    p[C.kprodmrnaikbe] = 1.052365e-01

    p[C.kshutnfkbin] = 2.072176e-01
    p[C.kshutnfkbout] = 5.537618e-04

    p[C.k0mrnaa20] = 1.118559e-03
    p[C.kdegmrnaa20] = 2.809754e-01
    p[C.kprodmrnaa20] = 1.788303e+00
    p[C.kpa20] = 2.445338e-01
    p[C.kdega20] = 6.217866e+00
    p[C.khillprodmrnaa20] = 1.243583e+00
    p[C.delayrnae] = 4.500000e+01

    p[C.X] = 0.0
    p[C.term] = 1.0

    return p
end


function initial_values()

    u0 = zeros(V.NUM)

    u0[V.B] = 1.0 # B
    u0[V.M] = 1.0 # M
    u0[V.C] = 1.0 # C
    u0[V.TAK1] = 1.0 # TAK1
    u0[V.IKK] = 1.0 # IKK
    u0[V.NFKBc] = 1.0 # NFKBc

    return u0
end