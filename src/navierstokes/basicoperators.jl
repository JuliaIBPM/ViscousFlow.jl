### Basic operators for any Navier-Stokes system ###
# Note that the input vector `w` is vorticity x cell spacing (i.e., it
# has units of velocity), and
# we expect the rhs of the equations to have units of dw/dt

export ns_rhs!

# Integrating factor -- rescale the time-step size
CartesianGrids.plan_intfact(t,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
       CartesianGrids.plan_intfact(t/(sys.Re*cellsize(sys)^2),w)

function ns_rhs!(dw::Nodes{Dual,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  dw .= 0.0
  _ns_rhs_convectivederivative!(dw,w,sys)
  _ns_rhs_double_layer!(dw,sys,t)
  return dw
end

function _ns_rhs_convectivederivative!(dw::Nodes{Dual,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY}
  Δx⁻¹ = 1/cellsize(sys)
  velocity!(sys.Vv,w,sys,0.0)
  _unscaled_convective_derivative!(sys.Vv,sys)
  sys.Sn .= 0.0
  curl!(sys.Sn,sys.Vv)
  sys.Sn .*= Δx⁻¹
  dw .-= sys.Sn
end

function _ns_rhs_double_layer!(dw::Nodes{Dual,NX,NY},
                              sys::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow},
                              t::Real) where {NX,NY,N,MT,FS}
  return dw
end

function _ns_rhs_double_layer!(dw::Nodes{Dual,NX,NY},
                              sys::NavierStokes{NX,NY,N,MT,FS,SD},
                              t::Real) where {NX,NY,N,MT,FS,SD}
  Δx⁻¹ = 1/cellsize(sys)
  fact = Δx⁻¹/sys.Re
  surface_velocity_jump!(sys.Δus,sys,t)
  sys.Vf .= 0.0
  sys.dlf(sys.Vf,sys.Δus)
  sys.Sn .= 0.0
  curl!(sys.Sn,sys.Vf)
  sys.Sn .*= fact
  dw .-= sys.Sn
end


#=
# RHS of Navier-Stokes (non-linear convective term)
function r₁(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY}) where {NX,NY,T}

 Ww = Edges(Dual,w) #sys.Ww
 Qq = Edges(Dual,w) #sys.Qq
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
=#
