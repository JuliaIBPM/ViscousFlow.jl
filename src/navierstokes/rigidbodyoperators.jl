#### Operators for a system with a body ####

# RHS of a stationary body with no surface velocity
function r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,N,T}
   ΔV = VectorData(sys.X̃)
   ΔV.u .-= sys.U∞[1]
   ΔV.v .-= sys.U∞[2]
   return ΔV
end

# Time-varying free stream
function r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true},U∞::RigidBodyTools.RigidBodyMotion) where {NX,NY,N,T}
   ΔV = VectorData(sys.X̃)
   _,ċ,_,_,_,_ = U∞(t)
   ΔV.u .-= real(ċ)
   ΔV.v .-= imag(ċ)
   return ΔV
end

# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ, B₂ = -ECL⁻¹
B₁ᵀ(f::VectorData{N},sys::NavierStokes{NX,NY,N,true}) where {NX,NY,N} = Curl()*(sys.Hmat*f)
B₂(w::Nodes{Dual,NX,NY,T},sys::NavierStokes{NX,NY,N,true}) where {NX,NY,T,N} = -(sys.Emat*(Curl()*(sys.L\w)))

# Constraint operators, using non-stored regularization and interpolation operators
B₁ᵀ(f::VectorData{N},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N} = Curl()*regop(sys.Fq,f)
B₂(w::Nodes{Dual,NX,NY,T},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,T,N} = -(regop(sys.Vb,Curl()*(sys.L\w)))

# Constraint operator constructors
# Constructor using stored operators
plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,T,N} =
                   (f -> B₁ᵀ(f,sys),w -> B₂(w,sys))

# Constructor using non-stored operators
function plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,T,N}
 regop = Regularize(sys.X̃,cellsize(sys);I0=origin(sys),issymmetric=true)

 return f -> B₁ᵀ(f,regop,sys),w -> B₂(w,regop,sys)
end




"""
    assign_velocity!(V::VectorData,X::VectorData,
                     xc::Real,yc::Real,α::Real,
                     mlist::Vector{RigidBodyMotion},t::Real)

Assign the components of rigid body velocity for every body (in inertial coordinate system)
at time `t` in the overall data structure `V`, using coordinates described by `X` (also in inertial
coordinate system), based on array of supplied motion `mlist` for each body.
"""
function RigidBodyTools.assign_velocity!(V::VectorData{N},X::VectorData{N},
                                           bl::BodyList,tlist::Vector{RigidTransform},
                                           mlist::Vector{RigidBodyMotion},t::Real) where {N}
    N == numpts(bl) || error("Inconsistent size of data structures")
    for i in 1:length(bl)
        ui = view(V.u,bl,i)
        vi = view(V.v,bl,i)
        xi = view(X.u,bl,i)
        yi = view(X.v,bl,i)
        Ti = tlist[i]
        assign_velocity!(ui,vi,xi,yi,Ti.trans[1],Ti.trans[2],Ti.α,mlist[i],t)
    end
end
