# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SimpleDiscontinuousGalerkin
# This script runs several minutes
trixi_include(joinpath(@__DIR__, "advection_linear_stability.jl"),
    N_elements = 20, tspan = (0.0, 3000.0),
    xrange = (-200, 10), yrange = (-200, 200), filename_extension = "_long_time",
    plot_solutions = false, plot_spectra = false)
