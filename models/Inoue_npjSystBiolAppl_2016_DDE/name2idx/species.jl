module V
const NAMES = [
  "B"
  "M"
  "BM"
  "C"
  "CB"
  "CM"
  "Cp"
  "CpB"
  "CpM"
  "CBM"
  "CpBM"
  #(* "TAK1 module*)
  "TAK1"
  "TAK1p"
  "TAK1C"
  "TAK1pC"
  #(* "IKK module*)
  "IKK"
  "IKKC"
  "IKKp"
  "IKKpC"
  "IKKppC"
  "IKKpp"
  "IKKi"
  #(* NFkB module*) c:cytoplasm n:nucleus
  "NFKBc"
  "NFKBn"
  "IKBac"
  "IKBbc"
  "IKBec"
  "IKBan"
  "IKBbn"
  "IKBen"
  "NFKBIKBac"
  "NFKBIKBbc"
  "NFKBIKBec"
  "NFKBIKBan"
  "NFKBIKBbn"
  "NFKBIKBen"
  "IKKIKBac"
  "IKKIKBbc"
  "IKKIKBec"
  "IKKNFKBIKBac"
  "IKKNFKBIKBbc"
  "IKKNFKBIKBec"
  "mRNAac"
  "mRNAbc"
  "mRNAec"
  "mRNAa20c"
  "A20c"
]

for (idx,name) in enumerate(NAMES)
  eval(Meta.parse("const $name = $idx"))
end

const NUM = length(NAMES)

end  # module