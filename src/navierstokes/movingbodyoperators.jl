#### Operators for a system with moving body ####

## Right-hand sides for rigid-body motion equations ##

# For a single rigid body
function r₁(u::Vector{Float64},t::Real,motion::RigidBodyMotion)
    _,ċ,_,_,α̇,_ = motion(t)
    return [real(ċ),imag(ċ),α̇]
end

# For systems of rigid bodies.
# The state vector is a concatenation of the 3 x 1 vectors for each body
function r₁(u::Vector{Float64},t::Real,ml::Vector{RigidBodyMotion})
    du = empty(u)
    z = zeros(Float64,3)
    for motion in ml
      append!(du,r₁(z,t,motion))
    end
    return du
    #du = similar(u)
    #cnt = 0
    #for ib = 1:length(ml)
    #    _,ċ,_,_,α̇,_ = ml[ib](t)
    #    du[cnt+1:cnt+3] = [real(ċ),imag(ċ),α̇]
    #    cnt += 3
    #end
    #return du
end

## Right-hand side for the combined Navier-Stokes and rigid-body motion equations ##

r₁(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY},
      motion::Union{RigidBodyMotion,Vector{RigidBodyMotion}}) where {NX,NY} =
     r₁(u[1],t,sys), r₁(u[2],t,motion)


## Constraint equations ##

function r₂(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false},
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

function plan_constraints(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N}

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
    return (f->B₁ᵀ(f,sys), f->zeros(Float64,size(u[2]))),
           (w->B₂(w,sys), u->Vector{Float64}())
  else
    return (f->B₁ᵀ(f,regop,sys), f->zeros(Float64,size(u[2]))),
           (w->B₂(w,regop,sys), u->Vector{Float64}())
  end



end
