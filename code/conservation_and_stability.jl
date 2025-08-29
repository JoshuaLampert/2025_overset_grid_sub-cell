# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SimpleDiscontinuousGalerkin
using LaTeXStrings
using Plots

EXAMPLES_DIR = joinpath(@__DIR__, "examples")
OUT = joinpath(@__DIR__, "figures")
isdir(OUT) || mkdir(OUT)
linewidth = 2

# Conservation and stability for advection, Maxwell's and Burgers' equations
plot_kwargs = (; linewidth=linewidth, title="", plot_title="", ylabel="", layout=(3, 1))
tspan = (0.0, 1.0)
trixi_include(joinpath(EXAMPLES_DIR, "subcell_advection_overset.jl"), tspan=tspan, interval=1, io=devnull)
p1 = plot(analysis_callback, exclude=(:mass, :entropy), label="advection", subplot=1; plot_kwargs...)
p2 = plot(analysis_callback, exclude=(:entropy, :entropy_timederivative), title="advection", legend=:outerright,
    label=L"\int \!u - \int u_0", subplot=1; plot_kwargs...)

trixi_include(joinpath(EXAMPLES_DIR, "subcell_maxwell_overset.jl"), tspan=tspan, interval=1, io=devnull)
plot!(p1, analysis_callback, exclude=(:electric_field, :magnetic_field, :entropy), label="Maxwell's", subplot=2; plot_kwargs...)
plot!(p2, analysis_callback, exclude=(:entropy, :entropy_timederivative), title="Maxwell's", legend=:outerright,
    linestyle=[:solid :dash], label=[L"\int \!E - \int E_0" L"\int \!B - \int B_0"], subplot=2; plot_kwargs...)

trixi_include(joinpath(EXAMPLES_DIR, "subcell_burgers_overset.jl"), tspan=tspan, interval=1, io=devnull, source_terms=nothing)
plot!(p1, analysis_callback, exclude=(:mass, :entropy), label="Burgers'", subplot=3; plot_kwargs...)
plot!(p2, analysis_callback, exclude=(:entropy, :entropy_timederivative), title="Burgers'", legend=:outerright,
    label=L"\int \!u - \int u_0", subplot=3; plot_kwargs...)

savefig(p1, joinpath(OUT, "subcell_overset_stability.pdf"))
savefig(p2, joinpath(OUT, "subcell_overset_conservation.pdf"))

# Compressible Euler equations (LLF loses conservation and stability, while HLL is conservative and stable)

n = 4
colors = palette(:default)[1:n]'
plot_kwargs_euler = (; linewidth=2, title="", plot_title="", ylabel="", layout=(2, 2))
linestyles_conservation = [:solid :dash :dot]
linestyle_stability = :dashdot
label_conservation = [L"\int\!\rho - \int\!\rho_0" L"\int\!\rho v - \int\!(\rho v)_0" L"\int\!\rho e - \int\!(\rho e)_0"]
label_stability = L"\int\!\partial S/\partial U \cdot U_t"

trixi_include(joinpath(EXAMPLES_DIR, "subcell_compressible_euler_overset.jl"), interval=1, io=devnull, source_terms=nothing)
p3 = plot(analysis_callback, exclude=(:entropy, :entropy_timederivative), label=label_conservation, subplot=1, legend=nothing,
    linestyles=linestyles_conservation, colors=colors[1:3]; plot_kwargs_euler...)
plot!(p3, analysis_callback, exclude=(:density, :momentum, :energy_total, :entropy), label=label_stability, subplot=2, legend=nothing,
    linestyles=linestyle_stability, color=colors[4]; plot_kwargs_euler...)

trixi_include(joinpath(EXAMPLES_DIR, "subcell_compressible_euler_overset.jl"), interval=1, io=devnull, source_terms=nothing, surface_flux=flux_lax_friedrichs)
plot!(analysis_callback, exclude=(:entropy, :entropy_timederivative), label=label_conservation, subplot=3, legend=nothing,
    linestyles=linestyles_conservation, colors=colors[1:3], legend_column=3; plot_kwargs_euler...)
plot!(p3, analysis_callback, exclude=(:density, :momentum, :energy_total, :entropy), label=label_stability, subplot=4, legend=nothing,
    linestyles=linestyle_stability, color=colors[4]; plot_kwargs_euler...)

plot!(p3, subplot = 3, legend = (-0.1, -0.6), bottom_margin = 15 * Plots.mm, legendfontsize = 8)
plot!(p3, subplot = 4, legend = (0.36, -0.6), legendfontsize = 8)
savefig(p3, joinpath(OUT, "subcell_overset_conservation_stability_euler.pdf"))
