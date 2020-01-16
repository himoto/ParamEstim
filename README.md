# Parameter estimation of ODE models describing biological processes

currently implimented for modeling early transcriptional regulation pathway ([Nakakuki *et al.*, ***Cell***, 2010](https://doi.org/10.1016/j.cell.2010.03.054))

![simulation_average](images/simulation_average.png)

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

## Model description
This mechanistic model describes the activation of immediate early genes such as cFos after epidermal growth factor (EGF) or heregulin (HRG) stimulation of the MAPK pathway. Phosphorylated cFos is a key transcription factor triggering downstream cascades of cell fate determination. The model can explain how the switch-like response of p-cFos emerges from the spatiotemporal dynamics. This mechanistic model comprises the explicit reaction kinetics of the signal transduction pathway, the transcriptional and the posttranslational feedback and feedforward loops. In the article, two different mechanistic models have been studied, the first one based on previously known interactions but failing to account for the experimental data and the second one including additional interactions which were discovered and confirmed by new experiments. The mechanistic model encoded here is the second one, the extended and at the time of creation most complete model of cell fate decision making in response to different doses of EGF or HRG stimulation.

## Dependencies
> - [Sundials](https://github.com/JuliaDiffEq/Sundials.jl)
> - [StatsBase](https://github.com/JuliaStats/StatsBase.jl)
> - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
> - [Seaborn](https://github.com/JuliaPy/Seaborn.jl)

## Usage
Parameter Estimation (*n*=1, 2, 3, · · ·)
```bash
$ mkdir logs
$ nohup julia optimize.jl n >> logs/n.log 2>&1 &
```
- If you want to continue from where you stopped in the last parameter search,
```bash
$ nohup julia optimize_continue.jl n >> logs/n.log 2>&1 &
```
---
Visualization of Simulation Results

    $ julia runSim.jl

```viz_type```:

- "average"
    : The average of simulation results with parameter sets in ```fitparam/```

- "best"
    : The best simulation result in ```fitparam/```, simulation with ```bestFitParam```

- "original"
    : Simulation with the default parameters and initial values defined in ```model/```

- "n(=1,2,...)"
    : Use the parameter set in ```fitparam/n/```

## Installation
    $ git clone https://github.com/himoto/ParamEstim.git


## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Kimura, S., Ono, I., Kita, H. & Kobayashi, S. An extension of UNDX based on guidelines for designing crossover operators: proposition and evaluation of ENDX. *Trans. Soc. Instrum. Control Eng.* **36**, 1162–1171 (2000). https://doi.org/10.9746/sicetr1965.36.1162

- Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance Independent Diversity Control for High Dimensional Function Optimization. *J. Japanese Soc. Artif. Intell.* **18**, 193–202 (2003). https://doi.org/10.1527/tjsai.18.193

- Kimura, S., Nakakuki, T., Kirita, S. & Okada, M. AGLSDC: A Genetic Local Search Suitable for Parallel Computation. *SICE J. Control. Meas. Syst. Integr.* **4**, 105–113 (2012). https://doi.org/10.9746/jcmsi.4.105

## License
[MIT](/LICENSE)
