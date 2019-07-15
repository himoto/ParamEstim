module V

const F_V = [
    "ppMEKc"
    "CREBn"
    "pCREBn"
    "ERKc"
    "ERKn"
    "pERKc"
    "pERKn"
    "ppERKc"
    "ppERKn"
    "Elk1n"
    "pElk1n"
    "cFOSc"
    "cFOSn"
    "pcFOSc"
    "pcFOSn"
    "DUSPc"
    "DUSPn"
    "pDUSPc"
    "pDUSPn"
    "DUSPn_ERKn"
    "DUSPn_pERKn"
    "DUSPn_ppERKn"
    "pDUSPn_ERKn"
    "pDUSPn_pERKn"
    "pDUSPn_ppERKn"
    "RSKc"
    "pRSKc"
    "pRSKn"
    "PrecfosmRNAn"
    "PreduspmRNAn"
    "cfosmRNAc"
    "duspmRNAc"
    "Fc"
    "Fn"
    "FmRNAc"
    "PreFmRNAn"
];

for (index,value) in enumerate(F_V)
    eval(Meta.parse("const $value=$index"));
end

const len_f_vars = length(F_V);

end  # module