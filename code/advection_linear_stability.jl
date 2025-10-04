# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SimpleDiscontinuousGalerkin
using LinearAlgebra
using Printf
using LaTeXStrings
using Plots

function stable(λ)
    tol = sqrt(eps(real(typeof(λ))))
    if real(λ) > tol
        return "unstable"
    else
        return "stable"
    end
end
function stable_markershape(λ)
    tol = sqrt(eps(real(typeof(λ))))
    if real(λ) > tol
        return :x
    else
        return :circle
    end
end

linewidth = 2
markersize = 2
EXAMPLES_DIR = joinpath(@__DIR__, "examples")
OUT = joinpath(@__DIR__, "figures")
isdir(OUT) || mkdir(OUT)

# We use an initial condition with a higher frequency to see the instability of classical
# overset grid method without sub-cell operators faster.
function initial_condition_convergence_test_higher_frequency(x, t, equations::LinearAdvectionEquation1D)
    x_trans = x - equations.advection_velocity * t
    return SVector(sinpi(4 * x_trans))
end
tspan_long = (0.0, 200.0)
a = -1.0
b = -0.1
c = 0.1
d = 1.0

p = 3

# We only use 9 elements for the left mesh to get the same number of DOFs compared to the case without sub-cell operators
mesh_left_subcell = Mesh(a, c, 9)
trixi_include(joinpath(EXAMPLES_DIR, "subcell_advection_overset.jl"), initial_condition=initial_condition_convergence_test_higher_frequency,
    a=a, b=b, c=c, d=d, p=p, mesh_left=mesh_left_subcell, tspan=tspan_long, interval=100, maxiters=1e9)
J = jacobian_fd(semi)
lamb_subcell = eigvals(J)
semi_subcell = semi
sol_subcell = sol
analysis_callback_subcell = analysis_callback

mesh_left = Mesh(a, c, 10)
Ds_left = [D_GLL for element in eachelement(mesh_left)]
solver_left = PerElementFDSBP(Ds_left,
    surface_integral=surface_integral,
    volume_integral=volume_integral)
trixi_include(joinpath(EXAMPLES_DIR, "subcell_advection_overset.jl"), initial_condition=initial_condition_convergence_test_higher_frequency,
    a=a, b=b, c=c, d=d, p=p, solver_left=solver_left, tspan=tspan_long, interval=100, maxiters=1e9)

plot(semi => sol, plot_title="", title="", ylabel="", label=["u" "v"], linewidth=linewidth, linestyle=[:dashdot :dashdotdot],
    color=[:purple :red])

plot!(semi_subcell => sol_subcell, plot_title="", title="", ylabel="", label=["u sub-cell" "v sub-cell"], linewidth=linewidth, linestyle=[:dash :dot],
    xticks=([a, -0.5, b, c, 0.5, d], ["a = $a", "-0.5", "b = $b", "c = $c", "0.5", "d = $d"]),
    color=[:blue :orange], legendfontsize=7)
mesh_fine = Mesh(a, d, 100)
x_fine = SimpleDiscontinuousGalerkin.left_element_boundary.(Ref(mesh_fine), 1:(nelements(mesh_fine)+1))
plot!(x_fine, only.(initial_condition.(x_fine, last(tspan), equations)), label="analytical solution", color=:black, linewidth=linewidth, linestyle=:solid)
# vline!(SimpleDiscontinuousGalerkin.left_element_boundary.(Ref(mesh_left), 1:(nelements(mesh_left)+1)),
#     color=:darkblue, alpha=0.5, label="", linestyle=:dash)
# Also plot the left mesh from the sub-cell run
vline!(SimpleDiscontinuousGalerkin.left_element_boundary.(Ref(mesh_left_subcell), 1:(nelements(mesh_left_subcell)+1)),
    color=:green, alpha=0.5, label="", linestyle=:dash)
vline!(SimpleDiscontinuousGalerkin.left_element_boundary.(Ref(mesh_right), 1:(nelements(mesh_right)+1)),
    color=:black, alpha=0.5, label="")

savefig(joinpath(OUT, "subcell_advection_overset.pdf"))

plot(analysis_callback, what=(:errors,), exclude=(:conservation_error,), title="", ylabel="Error", linewidth=linewidth,
    linestyles=[:solid :dash], label=[L"L^2" L"L^\infty"])
plot!(analysis_callback_subcell, what=(:errors,), exclude=(:conservation_error,), title="", ylabel="Error", linewidth=linewidth,
    linestyles=[:dashdot :dot], label=[L"$L^2$ sub-cell" L"$L^\infty$ sub-cell"])
savefig(joinpath(OUT, "subcell_advection_overset_errors.pdf"))

J = jacobian_fd(semi)
lamb = eigvals(J)

xrange = (-100, 5)
yrange = (-100, 100)
scatter(real.(lamb), imag.(lamb), group=stable.(lamb), xrange=xrange, yrange=yrange,
    title=@sprintf("max. real part: %.2e", maximum(real, lamb)),
    markersize=markersize, markershape=stable_markershape.(lamb),
    layout=2, subplot=1)
scatter!(real.(lamb_subcell), imag.(lamb_subcell), group=stable.(lamb_subcell), xrange=xrange, yrange=yrange,
    title=@sprintf("max. real part: %.2e", maximum(real, lamb_subcell)),
    markersize=markersize, markershape=stable_markershape.(lamb_subcell),
    layout=2, subplot=2)
savefig(joinpath(OUT, "subcell_advection_overset_spectra.pdf"))
