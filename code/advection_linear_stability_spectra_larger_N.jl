# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SimpleDiscontinuousGalerkin
using DoubleFloats: Double64
using Printf
using PrettyTables

# Create spectra for larger `N_elements`. We only use a small time span here since the spectra are constant in time due to the linearity of the problem.
# We use DoubleFloats here to get more accurate spectra since the eigenvalues can get positive because of the FD method used to compute the Jacobian.
tspan = (0.0, 1.0)
xranges = [(-75, 5), (-100, 5), (-200, 10), (-400, 20), (-800, 40)]
yranges = [(-75, 75), (-100, 100), (-200, 200), (-400, 400), (-800, 800)]
plot_solutions = false
plot_errors = false

table = []
for (i, N_elements) in enumerate([5, 10, 20, 40, 80])
    trixi_include(joinpath(@__DIR__, "advection_linear_stability.jl"),
        N_elements = N_elements, tspan = tspan, filename_extension = "_N$(N_elements)",
        xrange = xranges[i], yrange = yranges[i],
        RealT = Double64, plot_solutions = plot_solutions, plot_errors = plot_errors)
    lamb = @invokelatest (@__MODULE__).lamb
    lamb_subcell = @invokelatest (@__MODULE__).lamb_subcell
    push!(table, [string(N_elements), @sprintf("%.2e", maximum(real.(lamb))), @sprintf("%.2e", maximum(real.(lamb_subcell)))])
end
column_labels = ["N", "without sub-cell operator", "with sub-cell operator"]
style = LatexTableStyle(first_line_column_label=String[])
pretty_table(permutedims(hcat(table...)); column_labels, backend=:latex, table_format=latex_table_format__booktabs, style, alignment=:c)
