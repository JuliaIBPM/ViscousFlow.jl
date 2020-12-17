#### Operators for a system with moving body ####

export rigid_body_rhs!, rigid_body_rhs


## Right-hand sides for rigid-body motion equations ##

function rigid_body_rhs!(u::Vector{T},sys::NavierStokes,t::Real) where {T<:Real}
  length(u) == 3*(NDIM-1) || error("Wrong length for vector")
  u .= rigidbodyvelocity(sys.motions,t)
  return u
end

rigid_body_rhs(sys::NavierStokes,t::Real) where {T<:Real} = rigid_body_rhs!(zeros(Float64,3*(NDIM-1)),sys,t)

_construct_B₁ᵀ(sys::NavierStokes) = f->B₁ᵀ(f,sys)
_construct_B₂(sys::NavierStokes) = w->B₂(w,sys)
_construct_B₁ᵀ(u::Vector{T}) where {T<:Real} = zeros(u)
_construct_B₂(u::Vector{T}) where {T<:Real} = empty(u)


# For a single rigid body
#=
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

function r₂(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,MovingPoints},
                            motion::RigidBodyMotion) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.
  xc, yc, α = u[2]
  T = RigidBodyTools.RigidTransform((xc,yc),α)
  x, y = T(sys.points.u,sys.points.v)

  ΔV = VectorData(sys.points)
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

function plan_constraints(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,MovingPoints}) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.

  xc, yc, α = u[2]
  T = RigidBodyTools.RigidTransform((xc,yc),α)
  # should be able to save some time and memory allocation here...
  x, y = T(sys.points.u,sys.points.v)
  X = VectorData(x,y)
  regop = Regularize(X,cellsize(sys);I0=origin(g),weights=sys.areas.data,ddftype=ddftype)

  if sys._isstore
    Rf = RegularizationMatrix(regop,sys.Vb,sys.Ff)
    Ef = InterpolationMatrix(regop,sys.Ff,sys.Vb)
    sys.Rf = Rf
    sys.Ef = Ef
    return (f->B₁ᵀ(f,sys), f->zeros(Float64,size(u[2]))),
           (w->B₂(w,sys), u->Vector{Float64}())
  else
    return (f->B₁ᵀ(f,regop,sys), f->zeros(Float64,size(u[2]))),
           (w->B₂(w,regop,sys), u->Vector{Float64}())
  end

end
=#
