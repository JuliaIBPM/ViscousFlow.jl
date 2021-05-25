#### Operators for a system with a body and stationary points ####

export bc_constraint_rhs!

function bc_constraint_rhs!(us::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    # Here, need surface_velocity! - _velocity_single_layer!(u,sys,t) - _velocity_freestream!(u,sys,t)
    @unpack constraintop_cache = sys
    @unpack Vv, Δus = constraintop_cache
    Vv .= 0.0
    _velocity_single_layer!(Vv,sys,t)
    _velocity_freestream!(Vv,sys,t)
    us .= sys.Ef*Vv
    us .*= -1.0
    surface_velocity!(Δus,sys,t)
    us .+= Δus
    return us
end


# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ
function ns_op_constraint_force!(out::Nodes{Dual,NX,NY},τ::VectorData{N},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    @unpack constraintop_cache = sys
    @unpack Vv = constraintop_cache
    Vv .= sys.Rf*τ
    out .= 0.0
    curl!(out,Vv)
    return out
end

function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ,sys::NavierStokes{NX,NY,0}) where {NX,NY}
    return u
end

function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ::VectorData{N},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    u .= sys.Rf*τ
    return u
end


# B₂ = -ECL⁻¹
function bc_constraint_op!(out::VectorData{N},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    @unpack constraintop_cache = sys
    @unpack Vv = constraintop_cache
    Vv .= 0.0
    ViscousFlow._velocity_vorticity!(Vv,w,sys)
    out .= sys.Ef*Vv
    return out
end
