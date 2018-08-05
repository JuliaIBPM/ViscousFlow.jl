module Bodies

import Base:diff,length

export Body,RigidTransform,Ellipse,Plate

abstract type Body{N} end

"""
    RigidTransform(x::Tuple{Float64,Float64},α::Float64)

Construct a rigid-body transform operator, with rotation by angle `α` and
translation specified by `x`. The translation coordinates are specified in the
target coordinate system.

The resulting transform can be used as an operator on pairs of coordinate vectors,
`x` and `y`, or on bodies. For transformation of bodies, it only overwrites the
`x` and `y` fields of the body, but leaves the `x̃` and `ỹ` (body coordinates) intact.

The translation can be provided as either a tuple `(x,y)` or as a complex number.

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

RigidTransform(c::Complex128,α::Float64) = RigidTransform((real(c),imag(c)),α)

function Base.show(io::IO, T::RigidTransform)
    name = "Rigid-body transform"
    println(io, name)
    println(io, "  Translation: ($(T.trans[1]),$(T.trans[2]))")
    println(io, "  Rotation angle (rad): $(T.α)")
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
    diff(body::Body) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y differences of the faces on the perimeter of body `body`, whose centers
are at the current `x` and `y` coordinates (in inertial space) of the body.
"""
function diff(b::Body{N}) where {N}
  # need to modify this for thin flat plates

  ip1(i) = 1 + mod(i,N)
  im1(i) = 1 + mod(i-2,N)
  dxtmp = [0.5*(b.x[ip1(i)] - b.x[im1(i)]) for i = 1:N]
  dytmp = [0.5*(b.y[ip1(i)] - b.y[im1(i)]) for i = 1:N]
  return dxtmp,dytmp
end

"""
    dlength(body::Body) -> Vector{Float64}

Compute the lengths of the faces on the perimeter of body `body`, whose centers
are at the current `x` and `y` coordinates (in inertial space) of the body.
"""
dlength(b::Body{N}) where {N} = sqrt.(diff(b)[1].^2+diff(b)[2].^2)

"""
    normal(body::Body) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals (in inertial components) of the faces on the perimeter
of body `body`, whose centers are at the current `x` and `y` coordinates (in inertial space)
of the body.
"""
function normal(body::Body{N}) where {N}
  dx, dy = diff(body)
  ds = dlength(body)
  return -dy./ds, dx./ds
end

"""
    Ellipse(a,b,n)

Construct an elliptical body with semi-major axis `a` and semi-minor axis `b`,
with `n` points distributed on the body perimeter.

The constructor `Ellipse(a,n)` creates a circle of radius `a`.
"""
mutable struct Ellipse{N} <: Body{N}
  a :: Float64
  b :: Float64
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

end


function Ellipse(a::Float64,b::Float64,N::Int)
    x̃ = zeros(N)
    ỹ = zeros(N)
    θ = linspace(0,2π,N+1)
    @. x̃ = a*cos(θ[1:N])
    @. ỹ = b*sin(θ[1:N])


    Ellipse{N}(a,b,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ)
end

Ellipse(a::Float64,N::Int) = Ellipse(a,a,N)

function Base.show(io::IO, body::Ellipse{N}) where {N}
    println(io, "Elliptical body with $N points and semi-axes ($(body.a),$(body.b))")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

"""
    Plate(length,thick,n,[λ=1.0])

Construct a flat plate with length `length` and thickness `thick`,
with `n` points distributed on the body perimeter.

The optional parameter `λ` distributes the points differently. Values between `0.0`
and `1.0` are accepted.

The constructor `Plate(length,n,[λ=1.0])` creates a plate of zero thickness.
"""
mutable struct Plate{N} <: Body{N}
  len :: Float64
  thick :: Float64
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

end


function Plate(len::Float64,N::Int;λ::Float64=1.0)

    # set up points on plate
    #x = [[len*(-0.5 + 1.0*(i-1)/(N-1)),0.0] for i=1:N]

    Δϕ = π/(N-1)
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in linspace(π-Δϕ/2,Δϕ/2,N-1)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    x̃ = -0.5*len + Δϕ*cumsum([0.0; Jϕ])
    ỹ = zeros(x̃)

    Plate{N}(len,0.0,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ)

end

function Plate(len::Float64,thick::Float64,N::Int;λ::Float64=1.0)
    # input N is the number of panels on one side only

    # set up points on flat sides
    Δϕ = π/N
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in linspace(π-Δϕ/2,Δϕ/2,N)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    xtopface = -0.5*len + Δϕ*cumsum([0.0; Jϕ])
    xtop = 0.5*(xtopface[1:N] + xtopface[2:N+1])


    Δsₑ = Δϕ*Jϕ[1]
    Nₑ = 2*floor(Int,0.25*π*thick/Δsₑ)
    xedgeface = [0.5*len + 0.5*thick*cos(ϕ) for ϕ in linspace(π/2,-π/2,Nₑ+1)]
    yedgeface = [          0.5*thick*sin(ϕ) for ϕ in linspace(π/2,-π/2,Nₑ+1)]
    xedge = 0.5*(xedgeface[1:Nₑ]+xedgeface[2:Nₑ+1])
    yedge = 0.5*(yedgeface[1:Nₑ]+yedgeface[2:Nₑ+1])

    x̃ = Float64[]
    ỹ = Float64[]
    for xi in xtop
      push!(x̃,xi)
      push!(ỹ,0.5*thick)
    end
    for i = 1:Nₑ
      push!(x̃,xedge[i])
      push!(ỹ,yedge[i])
    end
    for xi in flipdim(xtop,1)
      push!(x̃,xi)
      push!(ỹ,-0.5*thick)
    end
    for i = Nₑ:-1:1
      push!(x̃,-xedge[i])
      push!(ỹ,yedge[i])
    end

    Plate{length(x̃)}(len,thick,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ)

end

function Base.show(io::IO, body::Plate{N}) where {N}
    println(io, "Plate with $N points and length $(body.len) and thickness $(body.thick)")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

#=
export Body

import Whirl
import Whirl.Grids

include("./RigidBodyMotions.jl")
using .RigidBodyMotions

struct BodyConfig

    "Position of body reference point in inertial space"
    xref::Array{Float64,1}

    "Rotation tensor of body"
    rot::Array{Float64,2}

end

function BodyConfig(xref::Vector{<:Real},angle::Float64)
    rot = [[cos(angle) -sin(angle)]; [sin(angle) cos(angle)]]

    BodyConfig(xref,rot)

end


mutable struct Body
    "number of Lagrange points on body"
    N::Int

    "coordinates of Lagrange points in body-fixed system"
    xtilde::Array{Array{Float64,1},1}

    "coordinates of Lagrange points in inertial system"
    x::Array{Array{Float64,1},1}

    "body motion function(s), which takes inputs t and position on body"
    motion::Vector{RigidBodyMotion}

    "body configuration"
    config::BodyConfig

end

# Set up blank body
Body() = Body(0,[])

function Body(N,xtilde)

    # set up array of inertial coordinates
    x = xtilde

    # default set the body motion function to 0 at every level
    motion = RigidBodyMotion[]
    push!(motion,RigidBodyMotion(0.0,0.0))

    # set configuration to origin and zero angle
    config = BodyConfig([0.0,0.0],0.0);

    Body(N,xtilde,x,motion,config)
end

function Body(N,xtilde,config::BodyConfig)

    b = Body(N,xtilde)
    update_body!(b,config)
    b

end

Body(N,xtilde,xref::Vector{Float64},angle::Float64) = Body(N,xtilde,BodyConfig(xref,angle))

function dims(body::Body)
    xmin = Inf*ones(Whirl.ndim)
    xmax = -Inf*ones(Whirl.ndim)
    for x in body.x
        xmin = [min(xmin[j],x[j]) for j = 1:Whirl.ndim]
        xmax = [max(xmax[j],x[j]) for j = 1:Whirl.ndim]
    end
    xmin, xmax
end

function update_body!(body::Body,config::BodyConfig)
    body.config = config
    body.x = transform(body.xtilde,config)
end

# Set the body motion function at level `l`
function set_velocity!(body::Body,m::RigidBodyMotion,l=1)
  if l > length(body.motion)
    push!(body.motion,m)
  else
    body.motion[l] = m
  end
end

# Evaluate some geometric details of a body




struct Identity{T, N, M}
end

function Base.show(io::IO, c::Identity{T,N,M}) where {T,N,M}
    print(io, "Identity operator for $M-dimensional vector field on $N body points")
end

function Identity(f::Array{T,2}) where {T}
  n,m = size(f)
  res = zeros(T,n,m)
  Identity{T,n,m}()
end

(c::Identity{T,N,M})(f) where {T,N,M} = f






function Base.show(io::IO, b::Body)
    println(io, "Body: number of points = $(b.N), "*
    		"reference point = ($(b.config.xref[1]),$(b.config.xref[2])), "*
		"rotation matrix = $(b.config.rot)")
    ds = map(x -> sqrt(x[1]^2 + x[2]^2),diff(b.x))
    println(io,"     max spacing between points = $(maximum(ds))")
    println(io,"     min spacing between points = $(minimum(ds))")
end

=#

end
