# This extends the package SimpleDiscontinuousGalerkin.jl for sub-cell operators in the context of
# overset grids. It provides a specialized surface integral for sub-cell operators.
# Note that we need to make an assumption on the wave speed direction. This has to be adjusted for negative wave speeds in the example script as follows:
#  - The example script has to use a sub-cell SBP operator on the `right_overlap_element` from the right mesh instead of
#    the `left_overlap_element` from the left mesh. Alternatively, one could also use a second sub-cell operator on the other mesh, but
#    it would integrate on the whole element, i.e., `WrappedSubcellOperator{:whole}(D)`, and not just on the overlap part of the element, i.e.,
#    `WrappedSubcellOperator{:left}(D)` or `WrappedSubcellOperator{:right}(D)`. This would essentially correspond to using one additional element.
# - The convex combination parameter `beta` in `SurfaceIntegralStrongFormSubcell` has to be set to `0.0` for left-traveling waves and to `1.0` for right-traveling waves.
# For problems with negative and positive wave speeds, it is even unclear what the mass is one wants to conserve because between the points b and c it is
# a convex combination of the integrals on the left and right mesh. It is unclear how to choose the convex combination parameter beta in general.

using SimpleDiscontinuousGalerkin
using SummationByPartsOperatorsExtra: SummationByPartsOperatorsExtra, AbstractNonperiodicDerivativeOperator, SubcellOperator, left_projection_right, right_projection_left, integrate_left, integrate_right

# We implement a wrapper type around a `SubcellOperator` to be able to add the information `PART` on whether to integrate on the left part on the right part of the element or on the whole element.
struct WrappedSubcellOperator{Part,T,Q,B,P,S} <: AbstractNonperiodicDerivativeOperator{T}
    D::SubcellOperator{T,Q,B,P,S}

    function WrappedSubcellOperator{Part}(D::SubcellOperator{T,Q,B,P,S}) where {Part,T,Q,B,P,S}
        new{Part,T,Q,B,P,S}(D)
    end
end

# This allows us to treat a `WrappedSubcellOperator` as a `SubcellOperator`. We forward all properties and all methods to the underlying `SubcellOperator`.
function Base.getproperty(D::WrappedSubcellOperator, name::Symbol)
    if name == :D
        return getfield(D, :D)
    else
        return getproperty(D.D, name)
    end
end

# It would be nicer to use some form of automatic method forwarding here, e.g., using MethodForwarding.jl, to avoid having to write out all the methods of `SubcellOperator` manually.
# However, there seem to be some issues with MethodForwarding.jl not forwarding all methods correctly. Instead, we manually forward the necessary methods here.
# using MethodForwarding: @forward
# @forward WrappedSubcellOperator{Part,T,Q,B,P,S} => SubcellOperator{T,Q,B,P,S} where {T,Q,B,P,S} [SummationByPartsOperatorsExtra]

Base.Matrix(D::WrappedSubcellOperator) = Matrix(D.D)
SummationByPartsOperatorsExtra.mass_matrix(D::WrappedSubcellOperator) = SummationByPartsOperatorsExtra.mass_matrix(D.D)
SummationByPartsOperatorsExtra.accuracy_order(D::WrappedSubcellOperator) = SummationByPartsOperatorsExtra.accuracy_order(D.D)
SummationByPartsOperatorsExtra.left_projection_right(D::WrappedSubcellOperator) = SummationByPartsOperatorsExtra.left_projection_right(D.D)
SummationByPartsOperatorsExtra.right_projection_left(D::WrappedSubcellOperator) = SummationByPartsOperatorsExtra.right_projection_left(D.D)
SummationByPartsOperatorsExtra.integrate(func, u, D::WrappedSubcellOperator) = SummationByPartsOperatorsExtra.integrate(func, u, D.D)

const SubcellOperatorOrWrapped = Union{SubcellOperator,WrappedSubcellOperator}

# For positive velocities/sub-cell operator on the left mesh we need `integrate_left` and
# for negative velocities/sub-cell operator on the right mesh we need `integrate_right`.
function SimpleDiscontinuousGalerkin.integrate_on_element(func, u, D::WrappedSubcellOperator{:left}, element, jacobian)
    return jacobian[element] * integrate_left(func, u, D.D)
end
function SimpleDiscontinuousGalerkin.integrate_on_element(func, u, D::WrappedSubcellOperator{:right}, element, jacobian)
    return jacobian[element] * integrate_right(func, u, D.D)
end
function SimpleDiscontinuousGalerkin.integrate_on_element(func, u, D::WrappedSubcellOperator{:whole}, element, jacobian)
    return jacobian[element] * integrate(func, u, D.D)
end

# For positive velocities/sub-cell operator on the left mesh use `left_projection_right`
# and `right_projection_left`, for negative velocities/sub-cell operator on the right mesh.
function SimpleDiscontinuousGalerkin.interpolation_operator(x, D::WrappedSubcellOperator{:left})
    return left_projection_right(D)
end
function SimpleDiscontinuousGalerkin.interpolation_operator(x, D::WrappedSubcellOperator{:right})
    return right_projection_left(D)
end

"""
    SurfaceIntegralStrongFormSubcell(beta, surface_integral::SurfaceIntegralStrongForm)

A specialized strong form surface integral for subcell operators. For a usual SBP operator,
this is equivalent to `SurfaceIntegralStrongForm`, but for subcell operators, it adds another
SAT term for the interface between the subcells.
"""
struct SurfaceIntegralStrongFormSubcell{T<:Real} <: SimpleDiscontinuousGalerkin.AbstractSurfaceIntegral
    surface_integral::SurfaceIntegralStrongForm
    beta::T
end

SurfaceIntegralStrongFormSubcell(beta) = SurfaceIntegralStrongFormSubcell(SurfaceIntegralStrongForm(), beta)

# This allows us to treat a `SurfaceIntegralStrongFormSubcell` as a `SurfaceIntegralStrongForm`.
function Base.getproperty(integral::SurfaceIntegralStrongFormSubcell, name::Symbol)
    if name in (:surface_integral, :beta)
        return getfield(integral, name)
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
        u_L = SimpleDiscontinuousGalerkin.get_node_vars(u, equations, 1, element) # u_a or v_b
        f_L = flux(u_L, equations)
        u_R = SimpleDiscontinuousGalerkin.get_node_vars(u, equations, nnodes(solver, element), element) # u_c or v_d
        f_R = flux(u_R, equations)
        surface_operator_left_ = SimpleDiscontinuousGalerkin.get_integral_operator(surface_operator_left, solver,
            element)
        surface_operator_right_ = SimpleDiscontinuousGalerkin.get_integral_operator(surface_operator_right, solver,
            element)
        D = SimpleDiscontinuousGalerkin.get_basis(solver, element)
        if D isa SubcellOperatorOrWrapped
            nvars = nvariables(equations)
            u_x_M_l = zeros(real(mesh), nvars) # u_{b_L} or v_{c_L}
            u_x_M_r = zeros(real(mesh), nvars) # u_{b_R} or v_{c_R}
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
            if D isa SubcellOperatorOrWrapped
                e_x_M_L = left_projection_right(D) # e_{b_L,u} or e_{c_L,v}
                e_x_M_R = right_projection_left(D) # e_{b_R,u} or e_{c_R,v}
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

# This only differs from the `calc_boundary_flux_left!` in SimpleDiscontinuousGalerkin.jl at the point c.
# For the sub-cell operator implementation, we need a convex combination of two SATs at the point c:
# -β * e_{c,u} * (f^{num}(u_c, v_{c_L}) - f(u_c)) - (1 - β) * e_{c,u} * (f^{num}(u_c, v_{c_R}) - f(u_c))
# = -e_{c,u} * (β * f^{num}(u_c, v_{c_L}) + (1 - β) * f^{num}(u_c, v_{c_R}) - f(u_c))
# Because we already compute
# -e_{c,u} * (surface_flux_values[v, 1, element] - f_L[v]) in the usual strong form surface integral (where f_L[v] is f(u_c)),
# see https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/blob/3fcfdabc7e54402e7ff8d96c0e2bb9cad76f2fc4/src/solvers/surface_integrals.jl#L111C35-L111C78,
# we only need to set surface_flux_values[v, 1, element] to β * f^{num}(u_c, v_{c_L}) + (1 - β) * f^{num}(u_c, v_{c_R}) instead of f^{num}(u_c, v_{c_L}).
function SimpleDiscontinuousGalerkin.calc_boundary_flux_left!(surface_flux_values_left, u, t, x_neg,
    equations, mesh, integral_left::SurfaceIntegralStrongFormSubcell, solver, cache)
    u_left, u_right = u
    solver_left, solver_right = solver
    mesh_left = mesh.mesh_left
    cache_left, _ = cache
    l_right = cache.l_right

    # Left boundary condition of left mesh (at a)
    e_left_L = SimpleDiscontinuousGalerkin.get_projection_operator(cache_left.e_left, solver_left, 1) # e_{a,u}
    u_ll = x_neg(u, SimpleDiscontinuousGalerkin.xmin(mesh), t, mesh, equations, solver, true, cache) # g_L
    u_rr = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_left, equations, e_left_L', :, 1) # u_a
    f = integral_left.surface_integral.surface_flux_boundary(u_ll, u_rr, equations)
    SimpleDiscontinuousGalerkin.set_node_vars!(surface_flux_values_left, f, equations, 1, 1)

    # Right boundary condition of left mesh (at c)
    D = SimpleDiscontinuousGalerkin.get_basis(solver_right, l_right)
    e_right_L = SimpleDiscontinuousGalerkin.get_projection_operator(cache_left.e_right, solver_left, nelements(mesh_left)) # e_{c,u}
    u_ll = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_left, equations, e_right_L', :, nelements(mesh_left)) # u_c
    # @assert D isa SubcellOperatorOrWrapped "The right overlap element of the left mesh should use a subcell operator for the surface integral to work correctly."
    if D isa SubcellOperatorOrWrapped
        beta = integral_left.beta
        e_M_left = left_projection_right(D) # e_{c_L,v}
        e_M_right = right_projection_left(D) # e_{c_R,v}
        u_rr_left = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_right, equations, e_M_left', :, l_right) # v_{c_L}
        u_rr_right = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_right, equations, e_M_right', :, l_right) # v_{c_R}
        f_left = integral_left.surface_integral.surface_flux_boundary(u_ll, u_rr_left, equations) # f^{num}(u_c, v_{c_L})
        f_right = integral_left.surface_integral.surface_flux_boundary(u_ll, u_rr_right, equations) # f^{num}(u_c, v_{c_R})
        f = beta * f_left + (1 - beta) * f_right # β * f^{num}(u_c, v_{c_L}) + (1 - β) * f^{num}(u_c, v_{c_R})
    else
        e_M_right = cache.e_M_right # e_{c,v}
        u_rr = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_right, equations, e_M_right', :, l_right) # v_c (v_{c_R})
        f = integral_left.surface_flux_boundary(u_ll, u_rr, equations) # f^{num}(u_c, v_c) (f^{num}(u_c, v_{c_R}))
    end
    SimpleDiscontinuousGalerkin.set_node_vars!(surface_flux_values_left, f, equations, 2, nelements(mesh_left))
    return nothing
end

# This only differs from the `calc_boundary_flux_left!` in SimpleDiscontinuousGalerkin.jl at the point b.
# See the comment above for the left boundary condition of the left mesh in `calc_boundary_flux_left!` for the explanation of the convex combination of the SATs.
function SimpleDiscontinuousGalerkin.calc_boundary_flux_right!(surface_flux_values_right, u, t, x_pos,
    equations, mesh, integral_right::SurfaceIntegralStrongFormSubcell, solver, cache)
    u_left, u_right = u
    solver_left, solver_right = solver
    mesh_right = mesh.mesh_right
    _, cache_right = cache
    l_left = cache.l_left

    # Left boundary condition of right mesh (at b)
    D = SimpleDiscontinuousGalerkin.get_basis(solver_left, l_left)
    e_left_R = SimpleDiscontinuousGalerkin.get_projection_operator(cache_right.e_left, solver_right, 1) # e_{b,v}
    u_rr = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_right, equations, e_left_R', :, 1) # v_b
    # @assert D isa SubcellOperatorOrWrapped "The left overlap element of the right mesh should use a subcell operator for the surface integral to work correctly."
    if D isa SubcellOperatorOrWrapped
        beta = integral_right.beta
        e_M_left = left_projection_right(D) # e_{b_L,u}
        e_M_right = right_projection_left(D) # e_{b_R,u}
        u_ll_left = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_left, equations, e_M_left', :, l_left) # u_{b_L}
        u_ll_right = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_left, equations, e_M_right', :, l_left) # u_{b_R}
        f_left = integral_right.surface_integral.surface_flux_boundary(u_ll_left, u_rr, equations) # f^{num}(u_{b_L}, v_b)
        f_right = integral_right.surface_integral.surface_flux_boundary(u_ll_right, u_rr, equations) # f^{num}(u_{b_R}, v_b)
        f = beta * f_left + (1 - beta) * f_right # β * f^{num}(u_{b_L}, v_b) + (1 - β) * f^{num}(u_{b_R}, v_b)
    else
        e_M_left = cache.e_M_left # e_{b,u}
        u_ll = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_left, equations, e_M_left', :, l_left) # u_b (u_{b_L})
        f = integral_right.surface_integral.surface_flux_boundary(u_ll, u_rr, equations) # f^{num}(u_b, v_b) (f^{num}(u_{b_L}, v_b))
    end
    SimpleDiscontinuousGalerkin.set_node_vars!(surface_flux_values_right, f, equations, 1, 1)

    # Right boundary condition of right mesh (at d)
    e_right_R = SimpleDiscontinuousGalerkin.get_projection_operator(cache_right.e_right, solver_right, nelements(mesh_right)) # e_{d,v}
    u_ll = SimpleDiscontinuousGalerkin.get_multiplied_node_vars(u_right, equations, e_right_R', :, nelements(mesh_right)) # v_d
    u_rr = x_pos(u, SimpleDiscontinuousGalerkin.xmax(mesh), t, mesh, equations, solver, false, cache) # g_R
    f = integral_right.surface_integral.surface_flux_boundary(u_ll, u_rr, equations)
    SimpleDiscontinuousGalerkin.set_node_vars!(surface_flux_values_right, f, equations, 2, nelements(mesh_right))
    return nothing
end

# This function defines how to combine the overlap integrals on the left and right mesh for the `AnalysisCallback`. For the sub-cell operator implementation, we use a convex combination of both integrals defined by the parameter `beta`.
function SimpleDiscontinuousGalerkin.overlap_integral_combination(surface_integral_left::SurfaceIntegralStrongFormSubcell, surface_integral_right,
    integral_overlap_left, integral_overlap_right)
    beta = surface_integral_left.beta
    integral_overlap = beta * integral_overlap_left + (1 - beta) * integral_overlap_right
    return integral_overlap
end
function SimpleDiscontinuousGalerkin.overlap_integral_combination(surface_integral_left, surface_integral_right::SurfaceIntegralStrongFormSubcell,
    integral_overlap_left, integral_overlap_right)
    beta = surface_integral_right.beta
    integral_overlap = beta * integral_overlap_left + (1 - beta) * integral_overlap_right
    return integral_overlap
end
function SimpleDiscontinuousGalerkin.overlap_integral_combination(surface_integral_left::SurfaceIntegralStrongFormSubcell, surface_integral_right::SurfaceIntegralStrongFormSubcell,
    integral_overlap_left, integral_overlap_right)
    beta_left = surface_integral_left.beta
    beta_right = surface_integral_right.beta
    @assert beta_left == beta_right "The beta parameters of the left and right surface integral should be the same for the overlap integral combination to work correctly."
    integral_overlap = beta_left * integral_overlap_left + (1 - beta_left) * integral_overlap_right
    return integral_overlap
end

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralStrongFormSubcell)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralStrongFormSubcell(", integral.surface_integral, ", beta = ", integral.beta, ")")
    end
end
