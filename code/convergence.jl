# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SimpleDiscontinuousGalerkin
import SummationByPartsOperatorsExtra
using PrettyTables
using Printf
using LaTeXStrings

EXAMPLES_DIR = joinpath(@__DIR__, "examples")

function create_eoc_table(filename)
    tables_str = []

    Ns = [10, 20, 50, 80]
    polydegs = [3, 4]


    for p in polydegs
        _, errorsmatrix = convergence_test(joinpath(EXAMPLES_DIR, filename), Ns, p=p, abstol=1e-14, reltol=1e-14)
        eocs = Dict(kind => log.(error[2:end, :] ./ error[1:(end-1), :]) ./
                            log.(Ns[1:(end-1)] ./ Ns[2:end])
                    for (kind, error) in errorsmatrix)
        table_str = []
        println(eocs)
        for (i, N) in enumerate(Ns)
            l2 = errorsmatrix[:l2][i, :]
            eoc_first_variable = i == 1 ? "---" : @sprintf("%.2f", eocs[:l2][i-1, 1])
            push!(table_str, [string(N), [@sprintf("%.2e", l2[v]) for v in eachvariable(equations)]..., eoc_first_variable])
        end
        push!(tables_str, table_str)
    end

    if equations isa CompressibleEulerEquations1D
        caption_equation_name = "compressible Euler equations"
        header_variables = [L"$\rho$", L"$\rho v$", L"$\rho e$"]
        label = "table:eocs_compressible_euler"
    else
        caption_equation_name = "linear advection equation"
        header_variables = ["error"]
        label = "table:eocs_advection"
    end
    column_labels = ["N", header_variables..., "EOC"]
    style = LatexTableStyle(first_line_column_label=String[])
    table_kwargs = (; column_labels, backend=:latex, table_format=latex_table_format__booktabs, style, alignment=:c)

    println("\\begin{table}[htb]")
    println("\\caption{\$L^2\$-errors and EOCs for the $caption_equation_name and polynomial degrees \$d=3\$ and \$d=4\$}")
    println("\\begin{subtable}{.5\\linewidth}")
    println("\\subcaption{\$d = 3\$}")
    println("\\centering")
    pretty_table(permutedims(hcat(tables_str[1]...)); table_kwargs...)
    println("\\end{subtable}%")
    println("\\begin{subtable}{.5\\linewidth}")
    println("\\subcaption{\$d = 4\$}")
    println("\\centering")
    pretty_table(permutedims(hcat(tables_str[2]...)); table_kwargs...)
    println("\\end{subtable}%")
    println("\\label{$label}")
    println("\\end{table}")
end

create_eoc_table("subcell_advection_overset.jl")
create_eoc_table("subcell_compressible_euler_overset.jl")
