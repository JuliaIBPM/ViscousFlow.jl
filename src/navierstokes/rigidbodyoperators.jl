#### Operators for a system with a body and stationary points ####

export bc_constraint_rhs!

function bc_constraint_rhs!(us::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    # Here, need surface_velocity! - _velocity_single_layer!(u,sys,t) - _velocity_freestream!(u,sys,t)
    sys.Vv .= 0.0
    _velocity_single_layer!(sys.Vv,sys,t)
    _velocity_freestream!(sys.Vv,sys,t)
    us .= sys.Ef*sys.Vv
    us .*= -1.0
    surface_velocity!(sys.Δus,sys,t)
    us .+= sys.Δus
    return us
end


# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ
function ns_op_constraint_force!(out::Nodes{Dual,NX,NY},τ::VectorData{N},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    sys.Vf .= sys.Rf*τ
    out .= 0.0
    curl!(out,sys.Vf)
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
    sys.Vv .= 0.0
    ViscousFlow._velocity_vorticity!(sys.Vv,w,sys)
    out .= sys.Ef*sys.Vv
    return out
end
