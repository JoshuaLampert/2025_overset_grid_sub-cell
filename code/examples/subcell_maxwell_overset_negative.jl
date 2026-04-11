using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK
using SummationByPartsOperatorsExtra: SubcellOperator, couple_subcell

include(joinpath("..", "surface_integral_subcell.jl"))

###############################################################################
# semidiscretization of the Maxwell equation

equations = MaxwellEquations1D(1.05)
function initial_condition_convergence_test_neg(x, t, equations::MaxwellEquations1D)
    c = equations.speed_of_light
    char_neg = x + c * t
    sin_char_neg = sinpi(2 * char_neg)

    E = -c * sin_char_neg
    B = sin_char_neg
    return SVector(E, B)
end
initial_condition = initial_condition_convergence_test_neg

a = -1.0
b = -0.1
c = 0.1
d = 1.0

N_elements = 10
mesh_left = Mesh(a, c, N_elements)
mesh_right = Mesh(b, d, N_elements)
mesh = OversetGridMesh(mesh_left, mesh_right)

x_L_ref = -1.0
x_R_ref = 1.0
p = 3
n_nodes = p + 1
D_GLL = legendre_derivative_operator(x_L_ref, x_R_ref, n_nodes)

# We only need a sub-cell operator for the right mesh because the waves travel to the left, but we could
# also additionally use `D_u_wrapped` at element `l_left` on the left mesh if we wanted to.
l_left = SimpleDiscontinuousGalerkin.left_overlap_element(mesh)
xl_L = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_left, l_left)
xl_R = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_left, l_left + 1)
linear_map(x, a, b, c, d) = c + (x - a) / (b - a) * (d - c)
b_mapped = linear_map(b, xl_L, xl_R, x_L_ref, x_R_ref)
D_left = legendre_derivative_operator(x_L_ref, b_mapped, n_nodes)
D_right = legendre_derivative_operator(b_mapped, x_R_ref, n_nodes)
D_u = couple_subcell(D_left, D_right, b_mapped)
D_u_wrapped = WrappedSubcellOperator{:whole}(D_u)

l_right = SimpleDiscontinuousGalerkin.right_overlap_element(mesh)
xr_L = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_right, l_right)
xr_R = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_right, l_right + 1)
c_mapped = linear_map(c, xr_L, xr_R, x_L_ref, x_R_ref)
D_left = legendre_derivative_operator(x_L_ref, c_mapped, n_nodes)
D_right = legendre_derivative_operator(c_mapped, x_R_ref, n_nodes)
D_v = couple_subcell(D_left, D_right, c_mapped)
D_v_wrapped = WrappedSubcellOperator{:right}(D_v)

surface_integral = SurfaceIntegralStrongForm(flux_godunov)
beta = 0.0 # beta = 0.0 corresponds to left-traveling waves
surface_integral_subcell = SurfaceIntegralStrongFormSubcell(surface_integral, beta)
volume_integral = VolumeIntegralStrongForm()

Ds_left = [D_GLL for element in eachelement(mesh_left)]
solver_left = PerElementFDSBP(Ds_left,
    surface_integral=surface_integral_subcell,
    volume_integral=volume_integral)
Ds_right = [element == l_right ? D_v_wrapped : D_GLL for element in eachelement(mesh_right)]
solver_right = PerElementFDSBP(Ds_right,
    surface_integral=surface_integral_subcell,
    volume_integral=volume_integral)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, (solver_left, solver_right);
    boundary_conditions=boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
tspan = (0.0, 1.0)
ode = SimpleDiscontinuousGalerkin.semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
io = stdout
analysis_callback = AnalysisCallback(semi; interval=100,
    extra_analysis_errors=(:conservation_error,), io=io)
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length=100)
sol = solve(ode, RDPK3SpFSAL49(), abstol=1e-8, reltol=1e-8,
    save_everystep=false, callback=callbacks, saveat=saveat)
