#module Streaming

import Base: *, +

using Interpolations
#using FastGaussQuadrature
using LinearAlgebra
using SpecialFunctions
using ForwardDiff
using DiffRules

import ForwardDiff:Dual,value,partials,derivative,extract_derivative

"""
    Params(ϵ,Re)

Set the parameters for the streaming solution. The problem is scaled so
that the radius of the cylinder is unity. Reynolds number is defined as Re = ΩR²/ν,
where Ω is the angular frequency of the oscillatory motion, e.g. sin(Ωt).
"""
struct Params
    ϵ :: Float64
    Re :: Float64
    γ² :: ComplexF64
    γ :: ComplexF64
    λ :: ComplexF64
    λ² :: ComplexF64
    H₀ :: ComplexF64
    C :: ComplexF64
  end

function Params(ϵ::Number,Re::Number)
        γ² = im*Re
        γ = exp(im*π/4)*√Re
        λ = √2*γ
        λ² = 2*γ²
        H₀ = hankelh1(0,γ)
        C = hankelh1(2,γ)/H₀
        Params(ϵ,Re,γ²,γ,λ,λ²,H₀,C)
    end

function Base.show(io::IO, p::Params) where {N}
        println(io, "Streaming flow parameters with Re = $(p.Re), ϵ = $(p.ϵ)")
end

"""
    ComplexFunc(f)

Provides a wrapper for a function expected to return complex values, for use in
dispatch in automatic differentiation with `ForwardDiff`.
"""
struct ComplexFunc{FT}
    fcn::FT
end

(f::ComplexFunc)(x) = f.fcn(x)

macro create_dual(f,nu,c,d,hankel)
    SF = :SpecialFunctions
    H = :($hankel)
    _f = Symbol("_",f)
    defs = quote
                $(_f)(x)  = $SF.$H($nu,$c*x)/$d
                $f = ComplexFunc($(_f))
           end
    return esc(defs)
end

macro extend_H(hankel)
    SF = :SpecialFunctions
    H = :($hankel)
    _, dH = DiffRules.diffrule(SF,H,:(ν),:(v))
    Hdef = quote
              function $SF.$H(ν,d::Complex{<:Dual{T}}) where {T}
                    dr, di = reim(d)
                    v = value(dr)+im*value(di)
                    return Dual{T}(real($SF.$H(ν,v)),real($dH)*partials(dr)-imag($dH)*partials(di)) +
                        im*Dual{T}(imag($SF.$H(ν,v)),imag($dH)*partials(dr)+real($dH)*partials(di))
              end
            end
    return esc(Hdef)
end

# Derivatives #

@inline function ForwardDiff.derivative(f::ComplexFunc{F}, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f.fcn, R))
    return ForwardDiff.extract_derivative(T, real(f.fcn(Dual{T}(x, one(x))))) +
        im*ForwardDiff.extract_derivative(T, imag(f.fcn(Dual{T}(x, one(x)))))
end

"""
    D²(f::ComplexFunc,K::Int) -> ComplexFunc

Return the radial part of the Laplacian at `r` of a function having cylindrical form

\$ f(r)e^{iKθ} \$
"""
@inline function D²(f::ComplexFunc,K::Int)
      rdf = ComplexFunc(r -> r*derivative(f,r))
      drdf = ComplexFunc(r -> derivative(rdf,r))
      return ComplexFunc(r -> drdf(r)/r - K^2*f(r)/r^2)
end

"""
    curl(f::ComplexFunc,K::Int) -> Tuple{ComplexFunc,ComplexFunc}

Return the radial part of the cylindrical components of the curl at `r` of a
function having cylindrical form

\$ f(r)e^{iKθ} \$
"""
@inline function curl(f::ComplexFunc,K::Int)
  return ComplexFunc(r -> K*f(r)/r),ComplexFunc(r -> -derivative(f,r))
end

# Integrals #
"""
    Integral(f,a<:Real,b<:Real[,fake_infinity=1e5],[,length=10000])

Set up an integral

\$ \\int_a^r f(\\tau)d\\tau \$

for a < r < b
If `b` is set to `Inf`, then this sets up the integral

\$ \\int_r^\\infty f(\\tau)d\\tau \$

where a < r < ∞
The optional argument `length` sets the size of the interpolation table, and
`fake_infinity` sets the maximum value at which to truncate the limit at
infinity.
"""
struct Integral{FT,V,IL}
    fcn::FT
    a::Float64
    b::Float64
    fake_infinity::Float64
    table::Interpolations.AbstractInterpolation
end

function Integral(fcn,a,b;fake_infinity=1e6,length=10000)
    V = typeof(fcn(1))
    v = zeros(V,length)
    τ = range(0,stop=1,length=length)
    if b == Inf
        # Integration variable τ = c/x - c/fake_infinity, where c is
        # such that τ goes from 1 to 0 when x goes from a to fake_infinity.
        lim_diff = 1/a-1/fake_infinity
        c = 1/lim_diff
        xi = c/(0+c/fake_infinity)
        fcn_old = fcn(xi)*xi^2/c
        for (i,τi) in enumerate(τ)
          if i > 1
            xi = c/(τi+c/fake_infinity)
            fcn_new = fcn(xi)*xi^2/c
            v[i] += v[i-1] + 0.5*step(τ)*(fcn_new+fcn_old)
            fcn_old = fcn_new
          end
        end
        table = LinearInterpolation(τ, v, extrapolation_bc = Flat())
        Integral{typeof(fcn),V,true}(fcn,a,b,fake_infinity,table)
    else
        # Integration variable τ = c(x - a), where c is such
        # that τ goes from 0 to 1 when x goes from a to b
        lim_diff = b-a
        c = 1/lim_diff
        fcn_old = fcn(a)
        for (i,τi) in enumerate(τ)
          if i > 1
            xi = a + τi*lim_diff
            fcn_new = fcn(xi)
            v[i] += v[i-1] + 0.5*step(τ)*lim_diff*(fcn_new+fcn_old)
            fcn_old = fcn_new
          end
        end
        table = LinearInterpolation(τ, v, extrapolation_bc = Flat())
        Integral{typeof(fcn),V,false}(fcn,a,b,fake_infinity,table)
    end
end

function Base.show(io::IO, I::Integral{FT,V,IL}) where {FT,V,IL}
    if IL
      println(io, "Integral from $(I.a) < r to $(I.b)")
    else
      println(io, "Integral from $(I.a) to r < $(I.b)")
    end
end

ComplexIntegral(fcn::ComplexFunc,a,b;ops...) = ComplexFunc(Integral(fcn,a,b;ops...))

ComplexIntegral(fcn::Function,a,b;ops...) = ComplexFunc(Integral(ComplexFunc(fcn),a,b;ops...))

function (I::Integral{FT,V,false})(r) where {FT,V<:Number}
  return r > I.a ? I.table((r-I.a)/(I.b-I.a)) : zero(V)
end

function (I::Integral{FT,V,true})(r) where {FT,V<:Number}
  return I.table((1/r-1/I.fake_infinity)/(1/I.a-1/I.fake_infinity))
end

#(I::Integral{FT,V,true} where {FT,V<:Number})(d::Dual{T}) where {T} = Dual{T}(I(value(d)),-I.fcn(value(d))*partials(d))

#(I::Integral{FT,V,false} where {FT,V<:Number})(d::Dual{T}) where {T} = Dual{T}(I(value(d)),I.fcn(value(d))*partials(d))

# Now define how they behave with Dual types. The following defines the automatic differentiation:
function (I::Integral{FT,V,true} where {FT<:ComplexFunc,V<:Number})(d::Dual{T}) where {T}
          return Dual{T}(real(I(value(d))),-real(I.fcn(value(d)))*partials(d)) +
       im*Dual{T}(imag(I(value(d))),-imag(I.fcn(value(d)))*partials(d))
end

(I::Integral{FT,V,false} where {FT<:ComplexFunc,V<:Number})(d::Dual{T}) where {T} =
         Dual{T}(real(I(value(d))),real(I.fcn(value(d)))*partials(d)) +
      im*Dual{T}(imag(I(value(d))),imag(I.fcn(value(d)))*partials(d))

# Set up functions #

@extend_H(hankelh1)
@extend_H(hankelh2)


@create_dual(X,0,p.γ,p.H₀,hankelh1)
@create_dual(Y,1,p.γ,p.H₀,hankelh1)
@create_dual(Z,2,p.γ,p.H₀,hankelh1)

@create_dual(H11,1,p.λ,1,hankelh1)
@create_dual(H12,1,p.λ,1,hankelh2)
@create_dual(H21,2,p.λ,1,hankelh1)
@create_dual(H22,2,p.λ,1,hankelh2)

abstract type Order end
abstract type First <: Order end
abstract type Second <: Order end

struct ComplexAmplitude{FT,K}
    f :: FT
end

ComplexAmplitude(f,K) = ComplexAmplitude{typeof(f),K == 1 ? First : Second}(f)

function (A::ComplexAmplitude{FT,Second})(x,y) where {FT}
    r2 = x^2+y^2
    r = sqrt(r2)
    return 2*real(A.f(r))*x*y/r2
end

#include("Utils.jl")
#using .Utils

# const nodes, weights = gausslegendre(100)
#
# function quadgauss(f::Function)
#     dot(weights,f(nodes))
# end
#
# struct Params
#   ϵ :: Float64
#   Re :: Float64
#   γ² :: Complex{Float64}
#   γ :: Complex{Float64}
#   λ :: Complex{Float64}
#   λ² :: Complex{Float64}
#   H₀ :: Complex{Float64}
#   X :: Function
#   Y :: Function
#   Z :: Function
#   C :: Complex{Float64}
# end
#
#
# function Params(ϵ::Float64,Re::Float64)
#   γ² = im*Re
#   γ = √γ²
#   λ = √2*γ
#   λ² = 2*γ²
#   H₀ = besselh(0,1,γ)
#   X(r) = besselh.(0,1,γ*r)./H₀
#   Y(r) = besselh.(1,1,γ*r)./H₀
#   Z(r) = besselh.(2,1,γ*r)./H₀
#   C = besselh.(2,1,γ)./H₀
#
#   Params(ϵ,Re,γ²,γ,λ,λ²,H₀,X,Y,Z,C)
# end
#
# struct Grid{N}
#   r :: Array{Float64,N}
#   Θ :: Array{Float64,N}
#   x :: Array{Float64,N}
#   y :: Array{Float64,N}
# end
#
# function Grid(x::Array{Float64,N},y::Array{Float64,N}) where {N}
#   r = similar(x)
#   Θ = similar(x)
#   @. r = sqrt(x^2+y^2)
#   @. Θ = atan2(y,x)
#   Grid(r,Θ,x,y)
# end
#
# Grid(x::Float64,y::Float64) = Grid([x],[y])
#
# function Base.show(io::IO, g::Grid{N}) where {N}
#     println(io, "$N-dimensional evaluation grid")
# end
#
# function D²(ψ::Vector{T},r::Vector{Float64},K::Int) where {T}
#   dψ = gradient(ψ,r)
#   gradient(r.*dψ,r)./r - K^2*ψ./r.^2
# end
#
# function curl(ψ::Vector{T},r::Vector{Float64},K::Int) where {T}
#   K*ψ./r,-gradient(ψ,r)
# end
#
#
#  struct ComplexAmplitude{N,K}
#   r :: Array{Float64,N}
#   ψ :: Array{Complex{Float64},N}
#   ω :: Array{Complex{Float64},N}
#   Ur :: Array{Complex{Float64},N}
#   UΘ :: Array{Complex{Float64},N}
# end
#
# function ComplexAmplitude(r::Vector{Float64},ψ::Vector{Complex{Float64}},K::Int)
#   ω = -D²(ψ,r,K)
#   Ur, UΘ = curl(ψ,r,K)
#   ComplexAmplitude{1,K}(r,ψ,ω,Ur,UΘ)
# end
#
# function Base.show(io::IO, s::ComplexAmplitude{N,K}) where {N,K}
#     #rmin = minimum(s.r)
#     #rmax = maximum(s.r)
#     println(io, "Amplitude of order $K streaming exact solution")
# end
#
# mutable struct Soln{T,N}
#   t :: T
#   ψ :: Array{T,N}
#   ω :: Array{T,N}
#   ur :: Array{T,N}
#   uθ :: Array{T,N}
# end
#
# function Soln(t::Float64,ψ0::Array{Complex{Float64},N},ω0::Array{Complex{Float64},N},
#           Ur0::Array{Complex{Float64},N},Uθ0::Array{Complex{Float64},N},K::Int) where {N}
#   ψ = real.(ψ0*exp(-im*K*t))
#   ω = real.(ω0*exp(-im*K*t))
#   ur = real.(Ur0*exp(-im*K*t))
#   uθ = real.(Uθ0*exp(-im*K*t))
#   Soln(t,ψ,ω,ur,uθ)
# end
#
# Soln(t::Float64,s₀::ComplexAmplitude{N,K}) where {N,K} = Soln(t,s₀.ψ,s₀.ω,s₀.ur,s₀.uΘ,K)
#
#
# function Base.show(io::IO, s::Soln{Float64,N}) where {N}
#     println(io, "Streaming exact solution at t = $(s.t)")
# end
#
# function Base.show(io::IO, s::Soln{Vector{Float64},N}) where {N}
#     println(io, "Streaming exact solution history from t = $(s.t[1]) to $(s.t[end])")
# end
#
# function Evaluate(t::Float64,g::Grid,s::ComplexAmplitude{N,K}) where {N,K}
#   @get s (r,ψ,ω,Ur,UΘ)
#
#   sin1 = sin.(K*g.Θ)
#   cos1 = cos.(K*g.Θ)
#
#   ψg = zeros(Complex{Float64},size(g.r))
#   ωg = zeros(Complex{Float64},size(g.r))
#   Urg = zeros(Complex{Float64},size(g.r))
#   UΘg = zeros(Complex{Float64},size(g.r))
#
#   ψitp = interpolate((r,),ψ,Gridded(Linear()))
#   ωitp = interpolate((r,),ω,Gridded(Linear()))
#   Uritp = interpolate((r,),Ur,Gridded(Linear()))
#   UΘitp = interpolate((r,),UΘ,Gridded(Linear()))
#
#   for (i,rd) in enumerate(g.r)
#     ψg[i] = ψitp[rd]*sin1[i]
#     ωg[i] = ωitp[rd]*sin1[i]
#     Urg[i] = Uritp[rd]*cos1[i]
#     UΘg[i] = UΘitp[rd]*sin1[i]
#   end
#   return Soln(t,ψg,ωg,Urg,UΘg,K)
#
# end
#
#
# function (a::Number * s::Soln)
#   snew = deepcopy(s)
#   snew.ψ *= a
#   snew.ω *= a
#   snew.ur *= a
#   snew.uθ *= a
#   snew
# end
#
# function (s1::Soln + s2::Soln)
#   snew = deepcopy(s1)
#   snew.ψ = s1.ψ + s2.ψ
#   snew.ω = s1.ω + s2.ω
#   snew.ur = s1.ur + s2.ur
#   snew.uθ = s1.uθ + s2.uθ
#   snew
# end
#
#
# Evaluate(g::Grid,s::ComplexAmplitude{N,K}) where {N,K} = Evaluate(0.0,g,s)
#
# function Evaluate(t::Float64,p::Params,g::Grid,s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {N}
#   s = p.ϵ*Evaluate(t,g,s₁)
#   s += p.ϵ^2*(Evaluate(g,s̄₂) + Evaluate(t,g,s₂))
#   s
# end
#
# function Evaluate1storder(t::Float64,p::Params,g::Grid,s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {N}
#   s = p.ϵ*Evaluate(t,g,s₁)
# #   s += p.ϵ^2*(Evaluate(g,s̄₂) + Evaluate(t,g,s₂))
#   s
# end
#
# function Evaluate2ndorder(t::Float64,p::Params,g::Grid,s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {N}
#   # s = p.ϵ*Evaluate(t,g,s₁)
#   s = p.ϵ^2*(Evaluate(g,s̄₂) + Evaluate(t,g,s₂))
#   s
# end
#
# function Evaluate(tr::Range{T},p::Params,g::Grid{N},s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {T,N}
#
# t = Float64[]
# ψ = [Float64[] for x in g.x]
# ω = [Float64[] for x in g.x]
# ur = [Float64[] for x in g.x]
# uθ = [Float64[] for x in g.x]
#
# for (i,ti) in enumerate(tr)
#   s = Evaluate(ti,p,g,s₁,s̄₂,s₂)
#   push!(t,ti)
#   for j = 1:length(g.x)
#     push!(ψ[j],s.ψ[j])
#     push!(ω[j],s.ω[j])
#     push!(ur[j],s.ur[j])
#     push!(uθ[j],s.uθ[j])
#   end
# end
#
# Soln(t,ψ,ω,ur,uθ)
#
# end
#
#
# function Evaluate1storder(tr::Range{T},p::Params,g::Grid{N},s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {T,N}
#
# t = Float64[]
# ψ = [Float64[] for x in g.x]
# ω = [Float64[] for x in g.x]
# ur = [Float64[] for x in g.x]
# uθ = [Float64[] for x in g.x]
#
# for (i,ti) in enumerate(tr)
#   s = Evaluate1storder(ti,p,g,s₁,s̄₂,s₂)
#   push!(t,ti)
#   for j = 1:length(g.x)
#     push!(ψ[j],s.ψ[j])
#     push!(ω[j],s.ω[j])
#     push!(ur[j],s.ur[j])
#     push!(uθ[j],s.uθ[j])
#   end
# end
#
# Soln(t,ψ,ω,ur,uθ)
#
# end
#
#
# function Evaluate2ndorder(tr::Range{T},p::Params,g::Grid{N},s₁::ComplexAmplitude{N,1},
#                             s̄₂::ComplexAmplitude{N,2},s₂::ComplexAmplitude{N,2}) where {T,N}
#
# t = Float64[]
# ψ = [Float64[] for x in g.x]
# ω = [Float64[] for x in g.x]
# ur = [Float64[] for x in g.x]
# uθ = [Float64[] for x in g.x]
#
# for (i,ti) in enumerate(tr)
#   s = Evaluate2ndorder(ti,p,g,s₁,s̄₂,s₂)
#   push!(t,ti)
#   for j = 1:length(g.x)
#     push!(ψ[j],s.ψ[j])
#     push!(ω[j],s.ω[j])
#     push!(ur[j],s.ur[j])
#     push!(uθ[j],s.uθ[j])
#   end
# end
#
# Soln(t,ψ,ω,ur,uθ)
#
# end
#
# function cartesian(s::Soln,g::Grid)
#   ux = s.ur.*cos.(g.Θ) - s.uθ.*sin.(g.Θ)
#   uy = s.ur.*sin.(g.Θ) + s.uθ.*cos.(g.Θ)
#   ux, uy
# end
#
# function cumint!(gint::Vector{T},g::Vector{T},x::Vector{Float64}) where {T}
#   gint[1] = 0.0
#   for i = 2:length(g)
#     gint[i] = gint[i-1] + 0.5*(g[i-1]+g[i])*(x[i]-x[i-1])
#   end
#   nothing
# end
#
# # integral of f₀(r) from r to ∞
# function frto∞(r,f₀)
#
#   tmp = similar(f₀.(r))
#   cumint!(tmp,f₀.(r[end:-1:1]).*r[end:-1:1].^2,1./r[end:-1:1])
#
#   # Compute the tail from rmax to ∞
#   rmax = maximum(r)
#   Irmaxto∞ = quadgauss() do x
#     if x == -1
#       return 0.0
#     else
#       t = 2*rmax./(x+1)
#       return 2*rmax.*f₀.(t)./(x+1).^2
#     end
#   end
#   tmp += Irmaxto∞
#
#   tmp[end:-1:1]
# end
#
# # integral of f₀(r)rᵅ from 1 to r
# function f1tor(r,f₀)
#   g = similar(f₀.(r))
#   cumint!(g,f₀.(r),r)
#   g
# end
#
#
#
# function FirstOrder(p::Params,r::Array{Float64,N}) where {N}
#   @get p (γ,X,Y,Z,C)
#   K = 1
#   ψ₀ = -(C./r -2Y.(r)./γ)
#   ψ̃₀ = ψ₀ - r
#   ComplexAmplitude(r,ψ₀,K)
# end
#
# function SecondOrderMean(p::Params,r::Array{Float64,N}) where {N}
#   @get p (Re,γ,γ²,X,Y,Z,C)
#
#   K = 2
#   f₀(r) = -0.5*γ²*Re*(0.5*(C*conj(X(r))-conj(C)*X(r))./r.^2 -
#                       0.5*conj(Z(r))+0.5*Z(r) +
#                       X(r).*conj(Z(r)) - conj(X(r)).*Z(r))
#
#
#   I⁻¹ = frto∞(r,r->f₀(r)/r)
#   I¹ = frto∞(r,r->f₀(r)*r)
#   I³ = f1tor(r,r->f₀(r)*r^3)
#   I⁵ = f1tor(r,r->f₀(r)*r^5)
#
#   ψ̃₀ = -r.^4/48.*I⁻¹ + r.^2/16.*I¹ +
#         I³/16 + I⁻¹[1]/16 - I¹[1]/8 +
#         1./r.^2.*(-I⁵/48-I⁻¹[1]/24+I¹[1]/16)
#
#   ψ₀ = ψ̃₀ - 0.5*im*(-C./r.^2+Z.(r))
#
#   ψsd = -0.5*imag((0-conj(X.(r))).*(C./r.^2-Z.(r)))
#   ψ₀ += ψsd
#
#   ComplexAmplitude(r,ψ₀,K)
# end
#
#
# function SecondOrder(p::Params,r::Array{Float64,N}) where {N}
#   @get p (Re,γ,γ²,λ,λ²,X,Y,Z,C)
#
#   K = 2
#   g₀(r) = 0.5*γ²*Re*(C*X(r)/r^2-Z(r));
#
#   H11(r) = besselh(1,1,λ*r);
#   H12(r) = besselh(1,2,λ*r);
#   H21(r) = besselh(2,1,λ*r);
#   H22(r) = besselh(2,2,λ*r);
#   Kλ(r) = H11(1)*H22(r) - H12(1)*H21(r);
#
#   IKgr = f1tor(r,r->r*Kλ(r)*g₀(r))
#   IH21gr = frto∞(r,r->r*H21(r)*g₀(r))
#   Igr⁻¹ = frto∞(r,r->g₀(r)/r)
#   Igr³ = f1tor(r,r->g₀(r)*r^3)
#
#   I¹ = 0.25*im*π/(λ²*H11(1))*IKgr.*H21.(r);
#   I² = 0.25*im*π/(λ²*H11(1))*IH21gr.*Kλ.(r);
#   I³ = 1/(λ²*λ*H11(1))*((H21.(r)-H21(1)./r.^2).*Igr⁻¹[1] + IH21gr[1]./r.^2);
#   I⁴ = -0.25/λ²*(Igr⁻¹.*r.^2-Igr⁻¹[1]./r.^2+Igr³./r.^2);
#   ψ̃₀ = I¹ + I² + I³ + I⁴;
#
#   ψ₀ = ψ̃₀ + 0.5*im*(-C./r.^2 + Z.(r));
#
#   ComplexAmplitude(r,ψ₀,K)
# end
#




#end
