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

#### Construct the first and second order solutions




#end
