# Parameter Estimation
Using Genetic Algorithm to Fit ODE Models to Data
***
![cfos](https://user-images.githubusercontent.com/31299606/50464653-81b02700-09d5-11e9-910a-e3e2dcbd4fdd.png)

 Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

## Usage
- Parameter Estimation
```julia
include("ParamEstim.jl")
using .ParamEstim
runGA()
```
- Visualization of Simulation Results
```julia
include("ParamEstim.jl")
using .ParamEstim
runSim()
```

## Installation
    $ git clone https://github.com/himoto/ParamEstim.git


## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Kimura, S., Ono, I., Kita, H. & Kobayashi, S. An extension of UNDX based on guidelines for designing crossover operators: proposition and evaluation of ENDX. *Trans. Soc. Instrum. Control Eng.* **36**, 1162–1171 (2000). https://doi.org/10.9746/sicetr1965.36.1162

- Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance Independent Diversity Control for High Dimensional Function Optimization. *J. Japanese Soc. Artif. Intell.* **18**, 193–202 (2003). https://doi.org/10.1527/tjsai.18.193

- Kimura, S., Nakakuki, T., Kirita, S. & Okada, M. AGLSDC: A Genetic Local Search Suitable for Parallel Computation. *SICE J. Control. Meas. Syst. Integr.* **4**, 105–113 (2012). https://doi.org/10.9746/jcmsi.4.105

## License
[MIT](/LICENSE)
