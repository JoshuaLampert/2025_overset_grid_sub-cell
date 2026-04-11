# Setup packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SummationByPartsOperatorsExtra

function print_all_matrices(D)
    println("D = ")
    display(Matrix(D))
    println("P = ")
    display(mass_matrix(D))
    println("B = ")
    display(mass_matrix_boundary(D))
    println("S = ")
    display(SummationByPartsOperatorsExtra.to_S(D))
    println()
end

ω_L = -1.0
ω_M = 0.0
ω_R = 1.0

# Example 5.4
n = 2
D_GLL_left = legendre_derivative_operator(ω_L, ω_M, n)
D_GLL_right = legendre_derivative_operator(ω_M, ω_R, n)
D_GLL_coupled = couple_subcell(D_GLL_left, D_GLL_right, ω_M)
println("coupled Gauss-Lobatto operators:")
print_all_matrices(D_GLL_coupled)

# Example 5.5
D_GR_left = polynomialbases_derivative_operator(GaussRadauLeft, ω_L, ω_M, n)
D_GR_right = polynomialbases_derivative_operator(GaussRadauRight, ω_M, ω_R, n)
D_GR_coupled = couple_subcell(D_GR_left, D_GR_right, ω_M)
println("coupled Gauss-Radau operators:")
print_all_matrices(D_GR_coupled)

# Example 5.6
N_L = 6
D_FD_left = derivative_operator(MattssonNordström2004(), 1, 2, ω_L, ω_M, N_L)
N_R = 8
D_FD_right = derivative_operator(MattssonNordström2004(), 1, 2, ω_M, ω_R, N_R)
D_FD_coupled = couple_subcell(D_FD_left, D_FD_right, ω_M)
println("coupled finite difference operators:")
print_all_matrices(D_FD_coupled)

# Example 5.7
D_GLL_FD_coupled = couple_subcell(D_GLL_left, D_FD_right, ω_M)
println("coupled Gauss-Lobatto finite difference operators:")
print_all_matrices(D_GLL_FD_coupled)
