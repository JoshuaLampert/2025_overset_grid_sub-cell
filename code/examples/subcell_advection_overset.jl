using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK
using SummationByPartsOperatorsExtra: SubcellOperator, couple_subcell

include(joinpath("..", "surface_integral_subcell.jl"))

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

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
l_left = SimpleDiscontinuousGalerkin.left_overlap_element(mesh)
xl_L = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_left, l_left)
xl_R = SimpleDiscontinuousGalerkin.left_element_boundary(mesh_left, l_left + 1)
linear_map(x, a, b, c, d) = c + (x - a) / (b - a) * (d - c)
b_mapped = linear_map(b, xl_L, xl_R, x_L_ref, x_R_ref)
D_left = legendre_derivative_operator(x_L_ref, b_mapped, n_nodes)
D_right = legendre_derivative_operator(b_mapped, x_R_ref, n_nodes)
D_u = couple_subcell(D_left, D_right, b_mapped)

surface_integral = SurfaceIntegralStrongForm(flux_godunov)
volume_integral = VolumeIntegralStrongForm()

Ds_left = [element == l_left ? D_u : D_GLL for element in eachelement(mesh_left)]
solver_left = PerElementFDSBP(Ds_left,
    surface_integral=SurfaceIntegralStrongFormSubcell(surface_integral),
    volume_integral=volume_integral)
Ds_right = [D_GLL for element in eachelement(mesh_right)]
solver_right = PerElementFDSBP(Ds_right,
    surface_integral=surface_integral,
    volume_integral=volume_integral)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, (solver_left, solver_right))

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 2.0
tspan = (0.0, 2.0)
ode = SimpleDiscontinuousGalerkin.semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
io = stdout
analysis_callback = AnalysisCallback(semi; interval=100,
    extra_analysis_errors=(:conservation_error,), io=io)
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length=100)
sol = solve(ode, RDPK3SpFSAL49(), abstol=1e-8, reltol=1e-8,
    save_everystep=false, callback=callbacks, saveat=saveat, maxiters=1e6)
