export Ellipse,Circle,Plate,NACA4


"""
    Ellipse(a,b,n) <: Body

Construct an elliptical body with semi-major axis `a` and semi-minor axis `b`,
with `n` points distributed on the body perimeter.
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
    θ = range(0,stop=2π,length=N+1)
    @. x̃ = a*cos(θ[1:N])
    @. ỹ = b*sin(θ[1:N])


    Ellipse{N}(a,b,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ)
end

"""
    Circle(a,n) <: Body

Construct a circular body with radius `a`
and with `n` points distributed on the body perimeter.
"""
Circle(a::Float64,N::Int) = Ellipse(a,a,N)

function Base.show(io::IO, body::Ellipse{N}) where {N}
    if body.a == body.b
      println(io, "Circular body with $N points and radius $(body.a)")
    else
      println(io, "Elliptical body with $N points and semi-axes ($(body.a),$(body.b))")
    end
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

"""
    Plate(length,thick,n,[λ=1.0]) <: Body

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
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in range(π-Δϕ/2,stop=Δϕ/2,length=N-1)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    x̃ = -0.5*len .+ Δϕ*cumsum([0.0; Jϕ])
    ỹ = zero(x̃)

    Plate{N}(len,0.0,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ)

end

function Plate(len::Float64,thick::Float64,N::Int;λ::Float64=1.0)
    # input N is the number of panels on one side only

    # set up points on flat sides
    Δϕ = π/N
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in range(π-Δϕ/2,stop=Δϕ/2,length=N)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    xtopface = -0.5*len .+ Δϕ*cumsum([0.0; Jϕ])
    xtop = 0.5*(xtopface[1:N] + xtopface[2:N+1])


    Δsₑ = Δϕ*Jϕ[1]
    Nₑ = 2*floor(Int,0.25*π*thick/Δsₑ)
    xedgeface = [0.5*len + 0.5*thick*cos(ϕ) for ϕ in range(π/2,stop=-π/2,length=Nₑ+1)]
    yedgeface = [          0.5*thick*sin(ϕ) for ϕ in range(π/2,stop=-π/2,length=Nₑ+1)]
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
    for xi in reverse(xtop,dims=1)
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

"""
    NACA4(cam,pos,thick[;np=20][,Xc=(0.0,0.0)][,len=1.0]) <: Body{N}

Generates points in the shape of a NACA 4-digit airfoil of chord length 1. The
relative camber is specified by `cam`, the position of
maximum camber (as fraction of chord) by `pos`, and the relative thickness
by `thick`.

The optional parameter `np` specifies the number of points on the upper
or lower surface. The optional parameter `Zc` specifies the mean position of
the vertices (which is set to the origin by default). The optional parameter
`len` specifies the chord length.

# Example

```jldoctest
julia> w = Bodies.NACA4(0.0,0.0,0.12);
```
"""
mutable struct NACA4{N} <: Body{N}
  len :: Float64
  camber :: Float64
  pos :: Float64
  thick :: Float64

  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

end


function NACA4(cam::Number,pos::Number,t::Number;np=20,Xc=(0.0,0.0),len=1.0)

# Here, cam is the fractional camber, pos is the fractional chordwise position
# of max camber, and t is the fractional thickness.

npan = 2*np-2

# Trailing edge bunching
an = 1.5
anp = an+1
x = zeros(np)

θ = zeros(size(x))
yc = zeros(size(x))

for j = 1:np
    frac = Float64((j-1)/(np-1))
    x[j] = 1 - anp*frac*(1-frac)^an-(1-frac)^anp;
    if x[j] < pos
        yc[j] = cam/pos^2*(2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam/pos*(1-x[j]/pos))
        end
    else
        yc[j] = cam/(1-pos)^2*((1-2*pos)+2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam*pos/(1-pos)^2*(1-x[j]/pos))
        end
    end
end

xu = zeros(size(x))
yu = xu
xl = xu
yl = yu

yt = t/0.20*(0.29690*sqrt.(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4)

xu = x-yt.*sin.(θ)
yu = yc+yt.*cos.(θ)

xl = x+yt.*sin.(θ)
yl = yc-yt.*cos.(θ)

rt = 1.1019*t^2;
θ0 = 0
if pos > 0
    θ0 = atan(2*cam/pos)
end
# Center of leading edge radius
xrc = rt*cos(θ0)
yrc = rt*sin(θ0)
θle = collect(0:π/50:2π)
xlec = xrc .+ rt*cos.(θle)
ylec = yrc .+ rt*sin.(θle)

# Assemble data
coords = [xu yu xl yl x yc]
cole = [xlec ylec]

# Close the trailing edge
xpanold = [0.5*(xl[np]+xu[np]); reverse(xl[2:np-1],dims=1); xu[1:np-1]]
ypanold = [0.5*(yl[np]+yu[np]); reverse(yl[2:np-1],dims=1); yu[1:np-1]]

xpan = zeros(npan)
ypan = zeros(npan)
for ipan = 1:npan
    if ipan < npan
        xpan1 = xpanold[ipan]
        ypan1 = ypanold[ipan]
        xpan2 = xpanold[ipan+1]
        ypan2 = ypanold[ipan+1]
    else
        xpan1 = xpanold[npan]
        ypan1 = ypanold[npan]
        xpan2 = xpanold[1]
        ypan2 = ypanold[1]
    end
    xpan[ipan] = 0.5*(xpan1+xpan2)
    ypan[ipan] = 0.5*(ypan1+ypan2)
end
w = ComplexF64[1;reverse(xpan,dims=1)+im*reverse(ypan,dims=1)]*len
w .-= mean(w)

x̃ = real.(w)
ỹ = imag.(w)


NACA4{length(x̃)}(len,cam,pos,t,Xc,0.0,x̃,ỹ,x̃,ỹ)

end

function Base.show(io::IO, body::NACA4{N}) where {N}
    println(io, "NACA 4-digit airfoil with $N points and length $(body.len) and thickness $(body.thick)")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end
