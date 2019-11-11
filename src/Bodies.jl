module Bodies

import Base:diff,length,push!,vec

using Statistics: mean

export Body,RigidTransform,midpoints,dlength,normal,dlengthmid,centraldiff,normalmid

abstract type Body{N} end

include("bodies/bodylist.jl")

"""
    RigidTransform(x::Tuple{Float64,Float64},α::Float64)

Construct a rigid-body transform operator, with rotation by angle `α` and
translation specified by `x`. The translation coordinates are specified in the
target coordinate system.

The resulting transform can be used as an operator on pairs of coordinate vectors,
`x` and `y`, or on bodies. For transformation of bodies, it only overwrites the
`x` and `y` fields of the body, but leaves the `x̃` and `ỹ` (body coordinates) intact.

The translation can be provided as either a tuple `(x,y)` or as a complex number.

# Constructors

- `RigidTransform((x,y),α)`
- `RigidTransform(u::Vector{Float64})`
- `RigidTransform(u::NTuple{3,Float64})`
- `RigidTransform.(u)` where `u` is a collection of vectors or tuples.

# Example

```jldoctest
julia> body = Bodies.Ellipse(0.5,0.1,100)
Elliptical body with 100 points and semi-axes (0.5,0.1)
   Current position: (0.0,0.0)
   Current angle (rad): 0.0

julia> T = RigidTransform((1.0,1.0),π/4)
Rigid-body transform
  Translation: (1.0,1.0)
  Rotation angle (rad): 0.7853981633974483

julia> T(body)
Elliptical body with 100 points and semi-axes (0.5,0.1)
   Current position: (1.0,1.0)
   Current angle (rad): 0.7853981633974483
```
"""
struct RigidTransform
   α   :: Float64
   rot :: Matrix{Float64}
   trans  :: Tuple{Float64,Float64}
end

function RigidTransform(x::Tuple{Float64,Float64},α::Float64)
    rot = [cos(α) -sin(α)
           sin(α) cos(α)]
    RigidTransform(α,rot,x)
end

RigidTransform(x::Union{Vector{Float64},NTuple{3,Float64}}) = RigidTransform((x[1],x[2]),x[3])

RigidTransform(c::ComplexF64,α::Float64) = RigidTransform((real(c),imag(c)),α)

#RigidTransform(x::Vector{Vector{Float64}}) =

function Base.show(io::IO, T::RigidTransform)
    name = "Rigid-body transform"
    println(io, name)
    println(io, "  Translation: ($(T.trans[1]),$(T.trans[2]))")
    println(io, "  Rotation angle (rad): $(T.α)")
end

"""
    vec(T::RigidTransform) -> Vector{Float64}

Returns a length-3 vector of the form [x,y,α] corresponding to the translation
and rotation specified by the given transform `T`.
"""
vec(T::RigidTransform) = [T.trans[1],T.trans[2],T.α]

"""
    vec(tl::Vector{RigidTransform}) -> NTuple{N,Vector{Float64}}

Returns a tuple of length-3 vectors of the form [x,y,α] corresponding to the translation
and rotation specified by the given by the list of transforms `tl`.
"""
function vec(tl::Union{Vector{RigidTransform},NTuple{N,RigidTransform}}) where {N}

  u = (vec(tl[1]),)
  for i = 2:length(tl)
    u = u...,vec(tl[i])
  end
  return u

end


function (T::RigidTransform)(x̃::Float64,ỹ::Float64)
    Xr = T.rot*[x̃,ỹ]
    return T.trans .+ (Xr[1],Xr[2])
end
function (T::RigidTransform)(x̃::AbstractVector{Float64},ỹ::AbstractVector{Float64})
    x = deepcopy(x̃)
    y = deepcopy(ỹ)
    for i = 1:length(x̃)
        x[i],y[i] = T(x̃[i],ỹ[i])
    end
    return x, y
end

function (T::RigidTransform)(b::Body{N}) where {N}
  b.x, b.y = T(b.x̃,b.ỹ)
  b.α = T.α
  b.cent = T.trans
  return b
end

# Evaluate some geometric details of a body
"""
    length(body::Body)

Return the number of points on the body perimeter
"""
length(::Body{N}) where {N} = N

"""
    diff(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y differences of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
"""
diff(b::Body{N}) where {N} = _diff(b.x,b.y)

function diff(bl::BodyList)
    dx = Float64[]
    dy = Float64[]
    for b in bl
        dxb, dyb = diff(b)
        append!(dx,dxb)
        append!(dy,dyb)
    end
    return dx, dy
end

function _diff(x::Vector{Float64},y::Vector{Float64})
  N = length(x)
  @assert N == length(y)

  ip1(i) = 1+mod(i,N)
  dxtmp = [x[ip1(i)] - x[i] for i = 1:N]
  dytmp = [y[ip1(i)] - y[i] for i = 1:N]

  return dxtmp, dytmp

end

"""
    midpoints(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y midpoints of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
"""
midpoints(b::Body{N}) where {N} = _midpoints(b.x,b.y)

function midpoints(bl::BodyList)
    xc = Float64[]
    yc = Float64[]
    for b in bl
        xcb, ycb = midpoints(b)
        append!(xc,xcb)
        append!(yc,ycb)
    end
    return xc, yc
end

function _midpoints(x::Vector{Float64},y::Vector{Float64})

  N = length(x)
  @assert N == length(y)

  ip1(i) = 1+mod(i,N)
  xc = 0.5*[x[ip1(i)] + x[i] for i = 1:N]
  yc = 0.5*[y[ip1(i)] + y[i] for i = 1:N]

  return xc, yc

end

"""
    centraldiff(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the circular central differences of coordinates on body `body` (or
on each body in list `body`).
"""
function centraldiff(b::Body{N}) where {N}
  ip1(i) = 1+mod(i,N)
  xc, yc = midpoints(b)
  xc .= circshift(xc,1)
  yc .= circshift(yc,1)

  return _diff(xc,yc)
end

function centraldiff(bl::BodyList)
    dx = Float64[]
    dy = Float64[]
    for b in bl
        dxb, dyb = centraldiff(b)
        append!(dx,dxb)
        append!(dy,dyb)
    end
    return dx, dy
end


"""
    dlength(body::Body/BodyList) -> Vector{Float64}

Compute the lengths of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example.
"""
function dlength(b::Union{Body,BodyList})
  dx, dy = diff(b)
  return sqrt.(dx.^2+dy.^2)
end


"""
    dlengthmid(body::Body/BodyList) -> Vector{Float64}

Compute the lengths of the faces formed between the face midpoints on the
perimeter of body `body`. The indexing of these midpoint faces is consistent
with that of the regular vertex points adjacent to both midpoints.
Midpoint face 2 corresponds to the face between midpoints 1 and 2, for example.
"""
function dlengthmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  return sqrt.(dx.^2+dy.^2)
end

"""
    normal(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals (in inertial components) of the faces on the perimeter
of body `body`, whose ends are at the current `x` and `y` coordinates (in inertial space)
of the body. Face 1 corresponds to the face between points 1 and 2, for example.
"""
function normal(b::Union{Body,BodyList})
  dx, dy = diff(b)
  ds = dlength(b)
  return -dy./ds, dx./ds
end

"""
    normalmid(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals (in inertial components) of the faces formed between
midpoints on the perimeter of body `body` (or each body in list `body`).
"""
function normalmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  ds = dlengthmid(b)
  return -dy./ds, dx./ds
end

include("bodies/shapes.jl")

end
