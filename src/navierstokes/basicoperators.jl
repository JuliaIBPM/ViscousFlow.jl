### Basic operators for any Navier-Stokes system ###

# Integrating factor -- rescale the time-step size
CartesianGrids.plan_intfact(Δt,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
       CartesianGrids.plan_intfact(Δt/(sys.Re*cellsize(sys)^2),w)

# RHS of Navier-Stokes (non-linear convective term)
function r₁(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY}) where {NX,NY,T}

 Ww = sys.Ww
 Qq = sys.Qq
 L = sys.L
 Δx⁻¹ = 1/cellsize(sys)

 grid_interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
 Qq.u .-= sys.U∞[1]
 Qq.v .-= sys.U∞[2]

 return rmul!(divergence(Qq∘grid_interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# RHS of Navier-Stokes (non-linear convective term)
function r₁(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY},U∞::RigidBodyTools.RigidBodyMotion) where {NX,NY,T}

 Ww = sys.Ww
 Qq = sys.Qq
 L = sys.L
 Δx⁻¹ = 1/cellsize(sys)

 grid_interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
 _,ċ,_,_,_,_ = U∞(t)
 Qq.u .-= real(ċ)
 Qq.v .-= imag(ċ)

 return rmul!(divergence(Qq∘grid_interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end
