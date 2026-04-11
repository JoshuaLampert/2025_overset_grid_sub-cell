# Numerical Experiments

This directory contains all code required to reproduce the numerical
experiments. First, you need to install Julia, e.g., by downloading
the binaries from the [download page](https://julialang.org/downloads/).
The numerical experiments were performed using Julia v1.12.5.

The code builds on the two Julia packages [SimpleDiscontinuousGalerkin.jl](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl)
(providing tools for 1D overset grid methods) and [SummationByPartsOperatorsExtra.jl](https://github.com/JoshuaLampert/SummationByPartsOperatorsExtra.jl)
(providing the implementation of sub-cell summation-by-parts operators).
The file `surface_integral_subcell.jl` provides glue code between these two libraries
implementing a surface integral tailored for sub-cell SBP operators.

The following list describes which script creates which figure(s) or tables
and the names of the resulting .pdf files:

* Figures 4(a), 4(b): `advection_linear_stability.jl` &rarr; `subcell_advection_overset.pdf`, `subcell_advection_overset_errors.pdf`
* Figure 5, 6(a), Table 2: `advection_linear_stability_spectra_larger_N.jl` &rarr; `subcell_advection_overset_spectra_N10.jl`, `subcell_advection_overset_spectra_N20.jl`
* Figures 6(b): `advection_linear_stability_long_time.jl` &rarr; `subcell_advection_overset_errors_long_time.pdf`
* Figures 7, 8, 9: `conservation_and_stability.jl` &rarr; `subcell_overset_conservation.pdf`, `subcell_overset_stability.pdf`, `subcell_overset_conservation_stability_euler.pdf`
* Tables 1, 3: `convergence.jl`

The resulting figures are then saved as .pdf files in a new directory `figures`
inside the folder of this `README.md`. The tables are printed to the screen as $\LaTeX$ code.

Additionally, the operators presented in Example 5.4 and Example 5.5 can be constructed and displayed using `examples.jl`.
Different polynomial degrees for which the operators are exact, can be used by changing the number of nodes
in each sub-cell `n`.

In order to execute a script, start Julia in this folder and execute

```julia
julia> include("file_name.jl")
```

in the Julia REPL. To execute the first script from the list above, e.g.,
execute

```julia
julia> include("advection_linear_stability.jl")
```
