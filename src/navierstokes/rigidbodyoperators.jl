#### Operators for a system with a body and stationary points ####

export constraint_rhs!

function constraint_rhs!(us::VectorData{N},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
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


# RHS of a stationary body with no surface velocity
#=
function r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,StaticPoints}) where {NX,NY,N,T}
   ΔV = VectorData(sys.points)
   ΔV.u .-= sys.U∞[1]
   ΔV.v .-= sys.U∞[2]
   return ΔV
end

# Time-varying free stream
function r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,StaticPoints},U∞::RigidBodyTools.RigidBodyMotion) where {NX,NY,N,T}
   ΔV = VectorData(sys.points)
   _,ċ,_,_,_,_ = U∞(t)
   ΔV.u .-= real(ċ)
   ΔV.v .-= imag(ċ)
   return ΔV
end
=#

# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ, B₂ = -ECL⁻¹
function B₁ᵀ(f::VectorData{N},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    sys.Vf .= sys.Rf*f
    sys.Sn .= 0.0
    curl!(sys.Sn,sys.Vf)
    return sys.Sn
end

function B₂(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    sys.Vv .= 0.0
    ViscousFlow._velocity_vorticity!(sys.Vv,w,sys)
    sys.Vb .= sys.Ef*sys.Vv
    return sys.Vb
end
#B₁ᵀ(f::VectorData{N},sys::NavierStokes{NX,NY,N,StaticPoints}) where {NX,NY,N} = Curl()*(sys.Rf*f)
#B₂(w::Nodes{Dual,NX,NY,T},sys::NavierStokes{NX,NY,N,StaticPoints}) where {NX,NY,T,N} = -(sys.Ef*(Curl()*(sys.L\w)))

# Constraint operators, using non-stored regularization and interpolation operators
#B₁ᵀ(f::VectorData{N},regop::Regularize,sys::NavierStokes{NX,NY,N,MovingPoints}) where {NX,NY,N} = Curl()*regop(sys.Ff,f)
#B₂(w::Nodes{Dual,NX,NY,T},regop::Regularize,sys::NavierStokes{NX,NY,N,MovingPoints}) where {NX,NY,T,N} = -(regop(sys.Vb,Curl()*(sys.L\w)))

# Constraint operator constructors
# Constructor using stored operators
plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,StaticPoints}) where {NX,NY,T,N} =
                   (f -> B₁ᵀ(f,sys),w -> B₂(w,sys))

# Constructor using non-stored operators
function plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,MovingPoints}) where {NX,NY,T,N}
 regop = Regularize(sys.points,cellsize(sys);I0=origin(sys),issymmetric=true)

 return f -> B₁ᵀ(f,regop,sys),w -> B₂(w,regop,sys)
end
