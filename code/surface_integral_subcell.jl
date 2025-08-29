# This extends the package SimpleDiscontinuousGalerkin.jl for subcell operators in the context of
# overset grids. It provides a specialized surface integral for subcell operators.
# Note that we need to make an assumption on the wave speed direction. It plays a role in three locations:
# 1. The specialization of `integrate_on_element` for subcell operators below (only relevant for the `AnalysisCallback`).
# 2. The choice of `left_projection_right` or `right_projection_left` in `interpolation_operator` below.
# 3. The choice whether the overlap region should be integrated on the left of on the right mesh. This is done
#    within SimpleDiscontinuousGalerkin.jl in the specialization of `PolynomialBasis.integrate` (only relevant
#    for the `AnalysisCallback`).

using SimpleDiscontinuousGalerkin
using SummationByPartsOperatorsExtra: SubcellOperator, left_projection_right, right_projection_left, integrate_left, integrate_right

# For positive velocities we need `integrate_left` and for negative velocities we need `integrate_right`.
function SimpleDiscontinuousGalerkin.integrate_on_element(func, u, D::SubcellOperator, element, jacobian)
    return jacobian[element] * integrate_left(func, u, D)
end

# For positive velocities use `left_projection_right`
# and `right_projection_left`, for negative velocities.
function SimpleDiscontinuousGalerkin.interpolation_operator(x, D::SubcellOperator)
    return left_projection_right(D)
end

"""
    SurfaceIntegralStrongFormSubcell(surface_integral::SurfaceIntegralStrongForm)

A specialized strong form surface integral for subcell operators. For a usual SBP operator,
this is equivalent to `SurfaceIntegralStrongForm`, but for subcell operators, it adds another
SAT term for the interface between the subcells.
"""
struct SurfaceIntegralStrongFormSubcell <: SimpleDiscontinuousGalerkin.AbstractSurfaceIntegral
    surface_integral::SurfaceIntegralStrongForm
end

SurfaceIntegralStrongFormSubcell() = SurfaceIntegralStrongFormSubcell(SurfaceIntegralStrongForm())

# This allows us to treat a `SurfaceIntegralStrongFormSubcell` as a `SurfaceIntegralStrongForm`.
function Base.getproperty(integral::SurfaceIntegralStrongFormSubcell, name::Symbol)
    if name == :surface_integral
        return getfield(integral, :surface_integral)
    else
        return getproperty(integral.surface_integral, name)
    end
end

function SimpleDiscontinuousGalerkin.create_cache(mesh, equations, solver, integral::SurfaceIntegralStrongFormSubcell)
    SimpleDiscontinuousGalerkin.create_cache(mesh, equations, solver, integral.surface_integral)
end

function SimpleDiscontinuousGalerkin.calc_surface_integral!(du, u, mesh, equations,
    surface_integral::SurfaceIntegralStrongFormSubcell, solver, cache)
    (; surface_operator_left, surface_operator_right, surface_flux_values) = cache
    for element in eachelement(mesh)
        u_L = SimpleDiscontinuousGalerkin.get_node_vars(u, equations, 1, element)
        f_L = flux(u_L, equations)
        u_R = SimpleDiscontinuousGalerkin.get_node_vars(u, equations, nnodes(solver, element), element)
        f_R = flux(u_R, equations)
        surface_operator_left_ = SimpleDiscontinuousGalerkin.get_integral_operator(surface_operator_left, solver,
            element)
        surface_operator_right_ = SimpleDiscontinuousGalerkin.get_integral_operator(surface_operator_right, solver,
            element)
        D = SimpleDiscontinuousGalerkin.get_basis(solver, element)
        if D isa SubcellOperator
            nvars = nvariables(equations)
            u_x_M_l = zeros(real(mesh), nvars)
            u_x_M_r = zeros(real(mesh), nvars)
            for v in eachvariable(equations)
                u_x_M_l[v] = left_projection_right(D)' * u[v, :, element]
                u_x_M_r[v] = right_projection_left(D)' * u[v, :, element]
            end
            f_L_subcell = flux(u_x_M_l, equations)
            f_R_subcell = flux(u_x_M_r, equations)
            surface_flux_value_subcell = surface_integral.surface_integral.surface_flux(u_x_M_l, u_x_M_r, equations)
        end
        for v in eachvariable(equations)
            du_update = surface_operator_left_ *
                        (surface_flux_values[v, 1, element] - f_L[v]) -
                        surface_operator_right_ *
                        (surface_flux_values[v, 2, element] - f_R[v])
            for node in eachnode(solver, element)
                du[v, node, element] += du_update[node]
            end
            # Additional SATs for subcell operator
            if D isa SubcellOperator
                e_x_M_L = left_projection_right(D)
                e_x_M_R = right_projection_left(D)
                inv_P = inv(mass_matrix(D))
                du_update = -inv_P * e_x_M_L * (surface_flux_value_subcell[v] - f_L_subcell[v]) +
                            inv_P * e_x_M_R * (surface_flux_value_subcell[v] - f_R_subcell[v])
                for node in eachnode(solver, element)
                    du[v, node, element] += du_update[node]
                end
            end
        end
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralStrongFormSubcell)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralStrongFormSubcell(", integral.surface_integral, ")")
    end
end
