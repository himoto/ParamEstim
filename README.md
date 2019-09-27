# Parameter Estimation
Using Genetic Algorithm to Fit ODE Models to Data

## Requirements
- **[Julia 1.0+](https://julialang.org)**
    - [IJulia](https://github.com/JuliaLang/IJulia.jl)
    - [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl)
    - [StatsBase](https://github.com/JuliaStats/StatsBase.jl)
    - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
    - [Seaborn](https://github.com/JuliaPy/Seaborn.jl)

## Usage
- Parameter Estimation (runGA/runGA_*n*.ipynb, *n*=1, 2, 3, · · ·)
```julia
display("text/html", """<script charset="utf-8">
    IPython.notebook.kernel.execute(
        'current_ipynb = "'+IPython.notebook.notebook_name+'" '
    );
    </script>"""
)
```
```julia
include("../ParamEstim.jl")
using .ParamEstim
optimize();

#= If you want to continue from where you stopeed in the last parameter search,

optimize_continue();

=#
```
- Visualization of Simulation Results (runSim.ipynb)
```julia
include("ParamEstim.jl");
using .ParamEstim
```
```julia
#==============================================================================
    viz_type::String => "best", "average", "original" or int(1~n_fitparam)
    show_all::Bool
    stdev::Bool (Only when viz_type == "average")
==============================================================================#

visualizeResult(Sim,viz_type="average",show_all=false,stdev=true)
```
1. **viz_type="best", show_all=true, stdev=false**

![sim_best](https://user-images.githubusercontent.com/31299606/61284133-54a05600-a7f9-11e9-93ea-2e249bcffd16.png)

2. **viz_type="average", show_all=false, stdev=true**

![sim_average](https://user-images.githubusercontent.com/31299606/61284148-5c5ffa80-a7f9-11e9-8dec-3b1f018f649a.png)

- Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.
## Installation
    $ git clone https://github.com/himoto/ParamEstim.git


## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Kimura, S., Ono, I., Kita, H. & Kobayashi, S. An extension of UNDX based on guidelines for designing crossover operators: proposition and evaluation of ENDX. *Trans. Soc. Instrum. Control Eng.* **36**, 1162–1171 (2000). https://doi.org/10.9746/sicetr1965.36.1162

- Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance Independent Diversity Control for High Dimensional Function Optimization. *J. Japanese Soc. Artif. Intell.* **18**, 193–202 (2003). https://doi.org/10.1527/tjsai.18.193

- Kimura, S., Nakakuki, T., Kirita, S. & Okada, M. AGLSDC: A Genetic Local Search Suitable for Parallel Computation. *SICE J. Control. Meas. Syst. Integr.* **4**, 105–113 (2012). https://doi.org/10.9746/jcmsi.4.105

## License
[MIT](/LICENSE)
