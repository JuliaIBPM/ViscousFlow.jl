module Bodies

import Base:diff

export Body,RigidTransform,Ellipse

abstract type Body{N} end


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

function (T::RigidTransform)(x̃::Float64,ỹ::Float64)
    Xr = T.rot*[x̃,ỹ]
    return T.trans .+ (Xr[1],Xr[2])
end
function (T::RigidTransform)(x̃::Vector{Float64},ỹ::Vector{Float64})
    x = copy(x̃)
    y = copy(ỹ)
    for i = 1:length(x̃)
        x[i],y[i] = T(x̃[i],ỹ[i])
    end
    return x, y
end
function (T::RigidTransform)(b::Body{N}) where {N}
  b.x, b.y = T(b.x̃,b.ỹ)
  b.α = T.α
  b.cent = T.trans
end

# Evaluate some geometric details of a body
function diff(b::Body{N}) where {N}

  ip1(i) = 1 + mod(i,N)
  im1(i) = 1 + mod(i-2,N)
  dxtmp = [0.5*(b.x̃[ip1(i)] - b.x̃[im1(i)]) for i = 1:N]
  dytmp = [0.5*(b.ỹ[ip1(i)] - b.ỹ[im1(i)]) for i = 1:N]
  return dxtmp,dytmp
end

ds(b::Body{N}) where {N} = sqrt.(diff(b)[1].^2+diff(b)[2].^2)


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

transform(x::Array{Float64,1},config::BodyConfig) = config.xref + config.rot*x

function transform(x::Array{Array{Float64,1},1},config::BodyConfig)
    [transform(x[i],config) for i=1:length(x)]
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


function normal(body::Body)
  nx = -dx(body)[2]./ds(body)
  ny = dx(body)[1]./ds(body)
  return [[nx[i],ny[i]] for i in 1:body.N]
end

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



function Ellipse(N::Int,a,b)::Body

    # set up the points on the circle with radius `rad`
    x = [[a*cos(2*pi*(i-1)/N),b*sin(2*pi*(i-1)/N)] for i=1:N]

    # put it at the origin, with zero angle
    Body(N,x,[0.0,0.0],0.0)

end

function Ellipse(N::Int,a,b,xcent::Vector{<:Real},angle)::Body

    b = Ellipse(N,a,b)
    update_body!(b,BodyConfig(xcent,angle))
    b

end

Circle(N::Int,rad) = Ellipse(N,rad,rad)

Circle(N::Int,rad,xcent::Vector{<:Real},angle) = Ellipse(N,rad,rad,xcent,angle)

function Plate(N::Int,len,λ)::Body

    # set up points on plate
    x = [[len*(-0.5 + 1.0*(i-1)/(N-1)),0.0] for i=1:N]

    Δϕ = π/N
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in linspace(π-Δϕ/2,Δϕ/2,N)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    xtop = -0.5*len + Δϕ*cumsum([0.0; Jϕ])

    x = [[xi,0.0] for xi in xtop];

    # put it at the origin, with zero angle
    Body(N,x,[0.0,0.0],0.0)

end

function Plate(N::Int,len,λ,xcent::Vector{<:Real},angle::Float64)::Body

    b = Plate(N,len,λ)
    update_body!(b,BodyConfig(xcent,angle))
    b

end

function Plate(N::Int,len,t,λ)::Body
    # input N is the number of panels on one side only

    # set up points on flat sides
    Δϕ = π/N
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in linspace(π-Δϕ/2,Δϕ/2,N)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    xtopface = -0.5*len + Δϕ*cumsum([0.0; Jϕ])
    xtop = 0.5*(xtopface[1:N] + xtopface[2:N+1])


    Δsₑ = Δϕ*Jϕ[1]
    Nₑ = 2*floor(Int,0.25*π*t/Δsₑ)
    xedgeface = [0.5*len + 0.5*t*cos(ϕ) for ϕ in linspace(π/2,-π/2,Nₑ+1)]
    yedgeface = [          0.5*t*sin(ϕ) for ϕ in linspace(π/2,-π/2,Nₑ+1)]
    xedge = 0.5*(xedgeface[1:Nₑ]+xedgeface[2:Nₑ+1])
    yedge = 0.5*(yedgeface[1:Nₑ]+yedgeface[2:Nₑ+1])


    x = [
         [[xi,0.5*t] for xi in xtop];
         [[xedge[i],yedge[i]] for i in 1:Nₑ];
         [[xi,-0.5*t] for xi in xtop[end:-1:1]];
         [[-xedge[i],yedge[i]] for i in Nₑ:-1:1]
        ]


    # put it at the origin, with zero angle
    Body(length(x),x,[0.0,0.0],0.0)

end

function Plate(N::Int,len,t,λ,xcent::Vector{<:Real},angle::Float64)::Body

    b = Plate(N,len,t,λ)
    update_body!(b,BodyConfig(xcent,angle))
    b

end

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
