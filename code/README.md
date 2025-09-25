# Numerical Experiments

This directory contains all code required to reproduce the numerical
experiments. First, you need to install Julia, e.g., by downloading
the binaries from the [download page](https://julialang.org/downloads/).
The numerical experiments were performed using Julia v1.11.6.

The code builds on the two Julia packages [SimpleDiscontinuousGalerkin.jl](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl)
(providing tools for 1D overset grid methods) and [SummationByPartsOperatorsExtra.jl](https://github.com/JoshuaLampert/SummationByPartsOperatorsExtra.jl)
(providing the implementation of sub-cell summation-by-parts operators).
The file `surface_integral_subcell.jl` provides glue code between these two libraries
implementing a surface integral tailored for sub-cell SBP operators.

The following list describes which script creates which figure(s) or tables
and the names of the resulting .pdf files:

* Figures 4(a), 4(b), 5: `advection_linear_stability.jl` &rarr; `subcell_advection_overset.pdf`, `subcell_advection_overset_errors.pdf`, `subcell_advection_overset_spectra.pdf`
* Figures 6, 7, 8: `conservation_and_stability.jl` &rarr; `subcell_overset_conservation.pdf`, `subcell_overset_stability.pdf`, `subcell_overset_conservation_stability_euler.pdf`
* Tables 1, 2: `convergence.jl`

The resulting figures are then saved as .pdf files in a new directory `figures`
inside the folder of this `README.md`. The tables are printed to the screen as $\LaTeX$ code.

In order to execute a script, start Julia in this folder and execute

```julia
julia> include("file_name.jl")
```

in the Julia REPL. To execute the first script from the list above, e.g.,
execute

```julia
julia> include("advection_linear_stability.jl")
```
