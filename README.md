# Parameter Estimation
Using Genetic Algorithm to Fit ODE Models to Data

## Requirements
- **[Julia 1.0+](https://julialang.org)**
    - [Sundials](https://github.com/JuliaDiffEq/Sundials.jl)
    - [StatsBase](https://github.com/JuliaStats/StatsBase.jl)
    - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
    - [Seaborn](https://github.com/JuliaPy/Seaborn.jl)

## Usage
Parameter Estimation (*n*=1, 2, 3, · · ·)
```bash
$ nohup julia optimize.jl n >> out/n.log 2>&1 &
```
- If you want to continue from where you stopped in the last parameter search,
```bash
$ nohup julia optimize_continue.jl n >> out/n.log 2>&1 &
```
---
Visualization of Simulation Results
```bash
$ julia runSim.jl
```

```viz_type```:

- "average"
    : The average of simulation results with parameter sets in ```out/```

- "best"
    : The best simulation result in ```out/```, simulation with ```best_fit_param```

- "original"
    : Simulation with the default parameters and initial values defined in ```biomass/model/```

- "n(=1,2,...)"
    : Use the parameter set in ```out/n/```

```julia
simulateAll(viz_type="best", show_all=true, stdev=false)
```

![simulation_best](images/simulation_best.png)

```julia
simulateAll(viz_type="average", show_all=false, stdev=true)
```

![simulation_average](images/simulation_average.png)

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

## Installation
    $ git clone https://github.com/himoto/ParamEstim.git


## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Kimura, S., Ono, I., Kita, H. & Kobayashi, S. An extension of UNDX based on guidelines for designing crossover operators: proposition and evaluation of ENDX. *Trans. Soc. Instrum. Control Eng.* **36**, 1162–1171 (2000). https://doi.org/10.9746/sicetr1965.36.1162

- Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance Independent Diversity Control for High Dimensional Function Optimization. *J. Japanese Soc. Artif. Intell.* **18**, 193–202 (2003). https://doi.org/10.1527/tjsai.18.193

- Kimura, S., Nakakuki, T., Kirita, S. & Okada, M. AGLSDC: A Genetic Local Search Suitable for Parallel Computation. *SICE J. Control. Meas. Syst. Integr.* **4**, 105–113 (2012). https://doi.org/10.9746/jcmsi.4.105

## License
[MIT](/LICENSE)
