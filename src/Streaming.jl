module Streaming

using Interpolations
using FastGaussQuadrature

include("Utils.jl")
using .Utils

const nodes, weights = gausslegendre(100)

function quadgauss(f::Function)
    dot(weights,f(nodes))
end

struct Params
  Re :: Float64
  γ² :: Complex{Float64}
  γ :: Complex{Float64}
  λ :: Complex{Float64}
  λ² :: Complex{Float64}
  H₀ :: Complex{Float64}
  X :: Function
  Y :: Function
  Z :: Function
  C :: Complex{Float64}
end

function Params(Re)
  γ² = im*Re
  γ = √γ²
  λ = √2*γ
  λ² = 2*γ²
  H₀ = besselh(0,1,γ)
  X(r) = besselh.(0,1,γ*r)./H₀
  Y(r) = besselh.(1,1,γ*r)./H₀
  Z(r) = besselh.(2,1,γ*r)./H₀
  C = besselh.(2,1,γ)./H₀

  Params(Re,γ²,γ,λ,λ²,H₀,X,Y,Z,C)
end

struct Grid{N}
  r :: Array{Float64,N}
  Θ :: Array{Float64,N}
  x :: Array{Float64,N}
  y :: Array{Float64,N}
end

function Grid(x::Array{Float64,N},y::Array{Float64,N}) where {N}
  r = similar(x)
  Θ = similar(x)
  @. r = sqrt(x^2+y^2)
  @. Θ = atan2(y,x)
  Grid(r,Θ,x,y)
end

function Base.show(io::IO, g::Grid{N}) where {N}
    println(io, "$N-dimensional evaluation grid")
end

function D²(ψ::Vector{T},r::Vector{Float64},order::Int) where {T}
  dψ = gradient(ψ,r)
  gradient(r.*dψ,r)./r - order^2*ψ./r.^2
end

function curl(ψ::Vector{T},r::Vector{Float64},order::Int) where {T}
  order*ψ./r,-gradient(ψ,r)
end


 struct ComplexAmplitude{N,K}
  r :: Array{Float64,N}
  ψ :: Array{Complex{Float64},N}
  ω :: Array{Complex{Float64},N}
  Ur :: Array{Complex{Float64},N}
  UΘ :: Array{Complex{Float64},N}
end

function ComplexAmplitude(r::Vector{Float64},ψ::Vector{Complex{Float64}},K::Int)
  ω = -D²(ψ,r,K)
  Ur, UΘ = curl(ψ,r,K)
  ComplexAmplitude{1,K}(r,ψ,ω,Ur,UΘ)
end

function Base.show(io::IO, s::ComplexAmplitude{N,K}) where {N,K}
    #rmin = minimum(s.r)
    #rmax = maximum(s.r)
    println(io, "Amplitude of order $K streaming exact solution")
end

struct Soln{N,K}
  t :: Float64
  ψ :: Array{Float64,N}
  ω :: Array{Float64,N}
  ur :: Array{Float64,N}
  uΘ :: Array{Float64,N}
end

function Soln(t::Float64,ψ0::Array{Complex{Float64},N},ω0::Array{Complex{Float64},N},
          Ur0::Array{Complex{Float64},N},UΘ0::Array{Complex{Float64},N},K::Int) where {N}
  ψ = real.(ψ0*exp(-im*K*t))
  ω = real.(ω0*exp(-im*K*t))
  ur = real.(Ur0*exp(-im*K*t))
  uΘ = real.(UΘ0*exp(-im*K*t))
  Soln{N,K}(t,ψ,ω,ur,uΘ)
end

Soln(t::Float64,s₀::ComplexAmplitude{N,K}) where {N,K} = Soln(t,s₀.ψ,s₀.ω,s₀.ur,s₀.uΘ,K)


function Base.show(io::IO, s::Soln{N,K}) where {N,K}
    println(io, "Order $K streaming exact solution at t = $(s.t)")
end

function Evaluate(t::Float64,g::Grid,s::ComplexAmplitude{N,K}) where {N,K}
  @get s (r,ψ,ω,Ur,UΘ)

  sin1 = sin.(K*g.Θ)
  cos1 = cos.(K*g.Θ)

  ψg = zeros(Complex{Float64},size(g.r))
  ωg = zeros(Complex{Float64},size(g.r))
  Urg = zeros(Complex{Float64},size(g.r))
  UΘg = zeros(Complex{Float64},size(g.r))

  ψitp = interpolate((r,),ψ,Gridded(Linear()))
  ωitp = interpolate((r,),ω,Gridded(Linear()))
  Uritp = interpolate((r,),Ur,Gridded(Linear()))
  UΘitp = interpolate((r,),UΘ,Gridded(Linear()))

  for (i,rd) in enumerate(g.r)
    ψg[i] = ψitp[rd]*sin1[i]
    ωg[i] = ωitp[rd]*sin1[i]
    Urg[i] = Uritp[rd]*cos1[i]
    UΘg[i] = UΘitp[rd]*sin1[i]
  end
  return Soln(t,ψg,ωg,Urg,UΘg,K)

end

Evaluate(g::Grid,s::ComplexAmplitude{N,K}) where {N,K} = Evaluate(0.0,g,s)

function cumint!(gint::Vector{T},g::Vector{T},x::Vector{Float64}) where {T}
  gint[1] = 0.0
  for i = 2:length(g)
    gint[i] = gint[i-1] + 0.5*(g[i-1]+g[i])*(x[i]-x[i-1])
  end
  nothing
end

# integral of f₀(r) from r to ∞
function frto∞(r,f₀)

  tmp = similar(f₀.(r))
  cumint!(tmp,f₀.(r[end:-1:1]).*r[end:-1:1].^2,1./r[end:-1:1])

  # Compute the tail from rmax to ∞
  rmax = maximum(r)
  Irmaxto∞ = quadgauss() do x
    if x == -1
      return 0.0
    else
      t = 2*rmax./(x+1)
      return 2*rmax.*f₀.(t)./(x+1).^2
    end
  end
  tmp += Irmaxto∞

  tmp[end:-1:1]
end

# integral of f₀(r)rᵅ from 1 to r
function f1tor(r,f₀)
  g = similar(f₀.(r))
  cumint!(g,f₀.(r),r)
  g
end



function FirstOrder(p::Params,r::Array{Float64,N}) where {N}
  @get p (γ,X,Y,Z,C)
  K = 1
  ψ₀ = -(C./r -2Y.(r)./γ)
  ψ̃₀ = ψ₀ - r
  ComplexAmplitude(r,ψ₀,K)
end

function SecondOrderMean(p::Params,r::Array{Float64,N}) where {N}
  @get p (Re,γ,γ²,X,Y,Z,C)

  K = 2
  f₀(r) = -0.5*γ²*Re*(0.5*(C*conj(X(r))-conj(C)*X(r))./r.^2 -
                      0.5*conj(Z(r))+0.5*Z(r) +
                      X(r).*conj(Z(r)) - conj(X(r)).*Z(r))


  I⁻¹ = frto∞(r,r->f₀(r)/r)
  I¹ = frto∞(r,r->f₀(r)*r)
  I³ = f1tor(r,r->f₀(r)*r^3)
  I⁵ = f1tor(r,r->f₀(r)*r^5)

  ψ̃₀ = -r.^4/48.*I⁻¹ + r.^2/16.*I¹ +
        I³/16 + I⁻¹[1]/16 - I¹[1]/8 +
        1./r.^2.*(-I⁵/48-I⁻¹[1]/24+I¹[1]/16)

  ψ₀ = ψ̃₀ - 0.5*im*(-C./r.^2+Z.(r))

  ψsd = -0.5*imag((0-conj(X.(r))).*(C./r.^2-Z.(r)))
  ψ₀ += ψsd

  ComplexAmplitude(r,ψ₀,K)
end


function SecondOrder(p::Params,r::Array{Float64,N}) where {N}
  @get p (Re,γ,γ²,λ,λ²,X,Y,Z,C)

  K = 2
  g₀(r) = 0.5*γ²*Re*(C*X(r)/r^2-Z(r));

  H11(r) = besselh(1,1,λ*r);
  H12(r) = besselh(1,2,λ*r);
  H21(r) = besselh(2,1,λ*r);
  H22(r) = besselh(2,2,λ*r);
  Kλ(r) = H11(1)*H22(r) - H12(1)*H21(r);

  IKgr = f1tor(r,r->r*Kλ(r)*g₀(r))
  IH21gr = frto∞(r,r->r*H21(r)*g₀(r))
  Igr⁻¹ = frto∞(r,r->g₀(r)/r)
  Igr³ = f1tor(r,r->g₀(r)*r^3)

  I¹ = 0.25*im*π/(λ²*H11(1))*IKgr.*H21.(r);
  I² = 0.25*im*π/(λ²*H11(1))*IH21gr.*Kλ.(r);
  I³ = 1/(λ²*λ*H11(1))*((H21.(r)-H21(1)./r.^2).*Igr⁻¹[1] + IH21gr[1]./r.^2);
  I⁴ = -0.25/λ²*(Igr⁻¹.*r.^2-Igr⁻¹[1]./r.^2+Igr³./r.^2);
  ψ̃₀ = I¹ + I² + I³ + I⁴;

  ψ₀ = ψ̃₀ + 0.5*im*(-C./r.^2 + Z.(r));

  ComplexAmplitude(r,ψ₀,K)
end


end
