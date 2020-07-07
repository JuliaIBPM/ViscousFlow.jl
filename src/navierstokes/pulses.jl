### Utilities ###

struct PointForce{T}
  x :: Vector{Float64}
  f0 :: Union{Float64,Vector{Float64}}
  t0 :: Float64
  σt :: Float64
  pulse :: PulseField
end

function PointForce(u::GridData,x0::Vector{Float64},f0,t0,σt,sys::NavierStokes)
  forcefield = GeneratedField(deepcopy(u),_createfields(f0,x0,sys),sys.grid)
  pulse = PulseField(forcefield,t0,σt)
  return PointForce{typeof(u)}(x0,f0,t0,σt,pulse)
end

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

#=
struct PointForce{T}
    x :: Tuple{Float64,Float64}
    f0 :: Union{Float64,Tuple{Float64,Float64}}
    σx :: Union{Nothing,Float64}
    t0 :: Float64
    σt :: Float64
    fbuffer :: Union{ScalarData,VectorData}
    xbuffer :: Union{Nothing,T}
    ubuffer :: T
    regop :: Regularize
end

"""
    PointForce(u::Union{Nodes,Edges},x0::Tuple{Float64,Float64},f0,t0,σt,sys::NavierStokes; σx)

Constructor function that immerses a point force in the `u`-type data of system `sys`,
of strength `f0` to be applied at physical position `x0`. This point force is modulated
in space by Gaussian with radial standard deviation `σx`(optional) and in time by a Gaussian centered at time `t0`
with standard deviation `σt`. The data `u` should be of either `Nodes` or `Edges` type.
If `Nodes`, then `f0` should be a scalar; if `Edges`, then `f0` should be a tuple.

The resulting function is a function of time and generates a field on `u`-type data.
"""
# PointForce version with Gaussian distribution in space (optional) and time centered about
# x0 and t0 with standard deviation σx in space and σt in time
function PointForce(w₀::T,x0::Tuple{Float64,Float64},
                    f0::Union{Float64,Tuple{Float64,Float64}},
                    t0::Float64,σt::Float64,sys::NavierStokes;
                    σx::Union{Nothing, Float64}=nothing) where {T <: Union{Nodes,Edges}}
    Xf = VectorData([x0[1]],[x0[2]])
    Ff = (T <: Nodes) ? ScalarData(Xf) : VectorData(Xf)
    Ff .= f0
    regop = Regularize(Xf,CartesianGrids.cellsize(sys);I0=CartesianGrids.origin(sys),issymmetric=true)
    is_spatial = typeof(σx)<:Float64
    if is_spatial
        xbuffer = zero(w₀)
        SpatialGauss(xbuffer, w₀, x0, σx, sys)
        return PointForce{T}(x0,f0,σx,t0,σt,Ff,xbuffer,T(),regop)
    else
        return PointForce{T}(x0,f0,nothing,t0,σt,Ff,nothing,T(),regop)
    end
end

# We add an Hadamard product of spatial and temporal Gaussians
function (f::PointForce)(t)
     is_spatial = typeof(f.σx)<:Float64
     if is_spatial
         return f.xbuffer∘f.regop(f.ubuffer,rmul!(deepcopy(f.fbuffer),exp(-(t-f.t0)^2/f.σt^2)))
     else
         return f.regop(f.ubuffer,rmul!(deepcopy(f.fbuffer),exp(-(t-f.t0)^2/f.σt^2)))
     end
end

function Base.show(io::IO,f::PointForce{T}) where {T}
    println(io,"Transient point force applied on the $(T) field.")
    println(io,"   strength = $(f.f0)")
    println(io,"   location = $(f.x)")
    if typeof(f.σx) <:Float64
    println(io,"   radial spread = $(f.σx)")
    end
    println(io,"   central time = $(f.t0)")
    println(io,"   half-interval = $(f.σt)")
end

struct Gaussian
  σ :: Real
  x0 :: Real
  y0 :: Real
end

(g::Gaussian)(r::Real) = exp(-0.5*r^2*g.σ^(-2))

#Gaussian!(r::Real, σ::Real) = exp(-0.5*r^2*σ^(-2))

#Gaussian(σ::Real) = r -> Gaussian!(r, σ)



function SpatialGauss(buffer::T, w₀::T,x0::Tuple{Float64,Float64},σx::Float64,sys::NavierStokes) where {T <: Union{Nodes,Edges}}

if T<: Nodes
        xx, yy = coordinates(w₀,dx=sys.grid.Δx,I0=origin(sys))
        for (i, xi) in enumerate(xx)
            for (j, yj) in enumerate(yy)
                buffer[i,j] = Gaussian(σx)(norm([xi-x0[1]; yj-x0[2]]))
            end
        end
else
    xx, xy, yx, yy = coordinates(w₀,dx=sys.grid.Δx,I0=origin(sys))
    for (i, xi) in enumerate(xx)
        for (j, yj) in enumerate(xy)
            buffer.u[i,j] = Gaussian(σx)(norm([xi-x0[1]; yj-x0[2]]))
        end
    end

    for (i, xi) in enumerate(yx)
        for (j, yj) in enumerate(yy)
            buffer.v[i,j] = Gaussian(σx)(norm([xi-x0[1]; yj-x0[2]]))
        end
    end
end

end
=#
