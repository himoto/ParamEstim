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
options[observables_index("Nuclear_NFkB")]["xlim"] = (0, 360)
options[observables_index("Nuclear_NFkB")]["xticks"] = [0, 120, 240,  360]
options[observables_index("Nuclear_NFkB")]["xlabel"] = "Time (min)"
options[observables_index("Nuclear_NFkB")]["ylim"] = (0.0, 1.1)
options[observables_index("Nuclear_NFkB")]["yticks"] = [0.0, 0.5, 1.0]
options[observables_index("Nuclear_NFkB")]["ylabel"] = "NFκB activity"

function set_rcParams()
    rc("figure",figsize = (4,3))
    rc("font",family = "Arial")
    rc("font",size = 20)
    rc("axes",linewidth = 1.5)
    rc("xtick.major",width = 1.5)
    rc("ytick.major",width = 1.5)
    rc("lines",linewidth = 1.8)
    rc("lines",markersize = 8)
end

# ---
end # module