struct PointForce{T}
  x :: Vector{Float64}
  f0 :: Union{Float64,Vector{Float64}}
  t0 :: Float64
  σt :: Float64
  pulse :: ModulatedField
end

function PointForce(u::GridData,x0::Vector{Float64},f0,t0,σt,sys::NavierStokes)
  forcefield = GeneratedField(deepcopy(u),_createfields(f0,x0,sys),sys.grid)
  pulse = PulseField(forcefield,t0,σt)
  return PointForce{typeof(u)}(x0,f0,t0,σt,pulse)
end

PointForce(u::GridData,x0::Tuple{Float64,Float64},f0,t0,σt,sys::NavierStokes) =
        PointForce(u,[x0...],f0,t0,σt,sys)


function _createfields(f0::Vector{T},x0,sys) where {T<:Real}
  fieldvec = AbstractSpatialField[]
  for f in f0
    push!(fieldvec,SpatialGaussian(cellsize(sys),x0...,f*cellsize(sys)^2))
  end
  return fieldvec
end

function _createfields(f0::T,x0,sys) where {T<:Real}
  return SpatialGaussian(cellsize(sys),x0...,f0*cellsize(sys)^2)
end

(p::PointForce)(t) = p.pulse(t)

function Base.show(io::IO,f::PointForce{T}) where {T}
    println(io,"Transient point force applied on the $(T) field.")
    println(io,"   strength = $(f.f0)")
    println(io,"   location = $(f.x)")
    println(io,"   central time = $(f.t0)")
    println(io,"   half-interval = $(f.σt)")
end
