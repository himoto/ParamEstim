module Viz
include("./observable.jl")

using PyPlot

const cm = PyPlot.cm.get_cmap("tab20")
options = [
    Dict(
        "divided_by" => 1.0,
        "xlim" => (),
        "xticks" => [],
        "xlabel" => "",
        "ylim" => (),
        "yticks" => [],
        "ylabel" => replace(replace(observables[i], "__" => "\n"), "_" => " "),
        "cmap" => [cm.colors[j] for j in 1:20],
        "shape" => ["o", "v", "^", "<", ">", "8", "s", 
                    "p", "*", "h", "H", "D", "d", "P", "X"],
    ) for i in 1:length(observables)]
# ---
for (i,_) in enumerate(observables)
    options[i]["divided_by"] = 60  # sec. -> min.
    options[i]["xlim"] = (-5, 95)
    options[i]["xticks"] = [0, 30, 60, 90]
    options[i]["xlabel"] = "Time (min)"
    options[i]["ylim"] = (-0.1, 1.3)
    options[i]["yticks"] = [0.0, 0.3, 0.6, 0.9, 1.2]
    options[i]["cmap"] = ["mediumblue", "red"]
    options[i]["shape"] = ["D", "s"]
end
options[
    observables_index("Phosphorylated_MEKc")
]["ylabel"] = "Phosphorylated MEK\n(cytoplasm)"
options[
    observables_index("Phosphorylated_ERKc")
]["ylabel"] = "Phosphorylated ERK\n(cytoplasm)"
options[
    observables_index("Phosphorylated_RSKw")
]["ylabel"] = "Phosphorylated RSK\n(whole cell)"
options[
    observables_index("Phosphorylated_CREBw")
]["ylabel"] = "Phosphorylated CREB\n(whole cell)"
options[
    observables_index("dusp_mRNA")
]["ylabel"] = L"$\it{dusp}$"*" mRNA\nexpression"
options[
    observables_index("cfos_mRNA")
]["ylabel"] = L"$\it{c}$"*"-"*L"$\it{fos}$"*" mRNA\nexpression"
options[
    observables_index("cFos_Protein")
]["ylabel"] = "c-Fos Protein\nexpression"
options[
    observables_index("Phosphorylated_cFos")
]["ylabel"] = "Phosphorylated c-Fos\nProtein expression"
# ---

function set_rcParams()
    rc("figure",figsize = (4,3))
    rc("font",family = "Arial")
    rc("font",size = 20)
    rc("axes",linewidth = 1.5)
    rc("xtick.major",width = 1.5)
    rc("ytick.major",width = 1.5)
    rc("lines",linewidth = 1.8)
    rc("lines",markersize = 12)
end

end # module