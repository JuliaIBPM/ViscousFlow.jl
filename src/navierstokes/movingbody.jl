function ConstrainedSystems.r₁(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY},
                            motion::RigidBodyMotion) where {NX,NY}
    _,ċ,_,_,α̇,_ = motion(t)
    return ConstrainedSystems.r₁(u[1],t,sys), [real(ċ),imag(ċ),α̇]
end

function ConstrainedSystems.r₂(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false},
                            motion::RigidBodyMotion) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.
  xc, yc, α = u[2]
  T = RigidBodyTools.RigidTransform((xc,yc),α)
  x, y = T(sys.X̃.u,sys.X̃.v)

  ΔV = VectorData(sys.X̃)
  _,ċ,_,_,α̇,_ = motion(t)
  for i = 1:N
      Δz = (x[i]-xc)+im*(y[i]-yc)
      ċi = ċ + im*α̇*Δz
      ΔV.u[i] = real(ċi)
      ΔV.v[i] = imag(ċi)
  end
  ΔV.u .-= sys.U∞[1]
  ΔV.v .-= sys.U∞[2]
  return ΔV, Vector{Float64}()
end

function ConstrainedSystems.plan_constraints(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.

  xc, yc, α = u[2]
  T = RigidBodyTools.RigidTransform((xc,yc),α)
  # should be able to save some time and memory allocation here...
  x, y = T(sys.X̃.u,sys.X̃.v)
  X = VectorData(x,y)
  regop = Regularize(X,cellsize(sys);issymmetric=true,I0=origin(sys))
  if sys._isstore
    Hmat, Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())
    sys.Hmat = Hmat
    sys.Emat = Emat
    return (f->ConstrainedSystems.B₁ᵀ(f,sys), f->zeros(Float64,size(u[2]))),
           (w->ConstrainedSystems.B₂(w,sys), u->Vector{Float64}())
  else
    return (f->ConstrainedSystems.B₁ᵀ(f,regop,sys), f->zeros(Float64,size(u[2]))),
           (w->ConstrainedSystems.B₂(w,regop,sys), u->Vector{Float64}())
  end



end
