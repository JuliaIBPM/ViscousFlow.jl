abstract type DDFType end

abstract type AbstractDDF end

struct DDF{C <: DDFType,OVERDX} <: AbstractDDF end

"""
    DDF([ddftype=Fields.Roma],[dx=1.0])

Construct a discrete delta function operator. This is generally only needed
internally by the `Regularize` operator, so the user doesn't have much need
for accessing this directly. The default DDF is the `Roma` function, which
has a support of 3 grid cells. Other choices are the `Goza` operator, which is
a truncated Gaussian with 28 cells support, and the `Witchhat`, which has 2
cells support. The resulting operator is evaluated with one, two or three
coordinate arguments, producing, respectively, 1-d, 2-d, or 3-d smeared delta functions.
It can also be called with the usual Julia vectorized dot notation with arrays of
arguments.
The optional cell spacing argument `dx` rescales the coordinates by this spacing,
and the result is also rescaled by this spacing (raised to the number of dimensions).
This spacing argument defaults to 1.0.

```jldoctest
julia> ddf = DDF(ddftype=Fields.Roma)
Discrete delta function operator of type ViscousFlow.Fields.Roma, with spacing 1.0

julia> ddf(1)
0.16666666666666666

julia> ddf(-1)
0.16666666666666666

julia> ddf.([-1,0,1])
3-element Array{Float64,1}:
 0.16666666666666666
 0.6666666666666666
 0.16666666666666666
```
"""
function DDF(;dx::Real=1.0,ddftype=Roma)
  DDF{ddftype,1.0/dx}()
end

struct GradDDF{C <: Fields.DDFType,OVERDX,D} <: AbstractDDF end

function GradDDF(dir::Int;dx::Real=1.0,ddftype=Roma)
  GradDDF{ddftype,1.0/dx,dir}()
end

# This macro creates an abstract type in the DDFType system for the named ddftype
# It also creates functions to evaluate the ddf in 1-d,2-d or 3-d.
macro ddffunc(ddftype)
    fname = Symbol("ddf_",lowercase(string(ddftype)))
    dfname = Symbol("ddf_d",lowercase(string(ddftype)))
    return esc(quote
            abstract type $ddftype <: DDFType end
            @inline (::DDF{$ddftype,OVERDX})(x::T) where {OVERDX,T <: Real} =
                  $fname(abs(x)*OVERDX)*OVERDX
            @inline (::DDF{$ddftype,OVERDX})(x::T,y::T) where {OVERDX,T <: Real} =
                          $fname(abs(x)*OVERDX)*OVERDX*$fname(abs(y)*OVERDX)*OVERDX
            @inline (::DDF{$ddftype,OVERDX})(x::T,y::T,z::T) where {OVERDX,T <: Real} =
                          $fname(abs(x)*OVERDX)*OVERDX*$fname(abs(y)*OVERDX)*OVERDX*$fname(abs(z)*OVERDX)*OVERDX

            @inline (::GradDDF{$ddftype,OVERDX,1})(x::T) where {OVERDX,T <: Real} =
                   $dfname(x*OVERDX)*OVERDX*OVERDX

            @inline (::GradDDF{$ddftype,OVERDX,1})(x::T,y::T) where {OVERDX,T <: Real} =
                   $dfname(x*OVERDX)*OVERDX*OVERDX*$fname(abs(y)*OVERDX)*OVERDX
            @inline (::GradDDF{$ddftype,OVERDX,2})(x::T,y::T) where {OVERDX,T <: Real} =
                   $fname(abs(x)*OVERDX)*OVERDX*$dfname(y*OVERDX)*OVERDX*OVERDX

            @inline (::GradDDF{$ddftype,OVERDX,1})(x::T,y::T,z::T) where {OVERDX,T <: Real} =
                   $dfname(x*OVERDX)*OVERDX*OVERDX*$fname(abs(y)*OVERDX)*OVERDX*$fname(abs(z)*OVERDX)*OVERDX
            @inline (::GradDDF{$ddftype,OVERDX,2})(x::T,y::T,z::T) where {OVERDX,T <: Real} =
                   $fname(abs(x)*OVERDX)*OVERDX*$dfname(y*OVERDX)*OVERDX*OVERDX*$fname(abs(z)*OVERDX)*OVERDX
            @inline (::GradDDF{$ddftype,OVERDX,3})(x::T,y::T,z::T) where {OVERDX,T <: Real} =
                   $fname(abs(x)*OVERDX)*OVERDX*$fname(abs(y)*OVERDX)*OVERDX*$dfname(z*OVERDX)*OVERDX*OVERDX
            end)
end

#==== ROMA ====#

@ddffunc Roma

roma1(r) = (1+sqrt(-3r^2+1))/3
roma2(r) = (5-3r-sqrt(1-3(1-r)^2))/6

@inline ddf_roma(r::Real) = r > 1.5 ? 0.0 : r <= 0.5 ? roma1(r) : roma2(r)

droma1(r) = -r/sqrt(-3r^2+1)
droma2(r) = -0.5*(1+(1-r)/sqrt(1-3(1-r)^2))

@inline ddf_droma(r::Real) = abs(r) > 1.5 ? 0.0 : abs(r) <= 0.5 ? sign(r)*droma1(abs(r)) : sign(r)*droma2(abs(r))

#==== GOZA ====#

@ddffunc Goza

@inline ddf_goza(r::Real) = r > 14 ? 0.0 : exp(-π^2/36 * r^2)*sqrt(π/36)

@inline ddf_dgoza(r::Real) = abs(r) > 14 ? 0.0 : (-π^2/18)*r*exp(-π^2/36 * r^2)*sqrt(π/36)

#==== WITCHHAT ====#

@ddffunc Witchhat

@inline ddf_witchhat(r::Real) = r > 1 ? 0.0 : 1-r

@inline ddf_dwitchhat(r::Real) = abs(r) > 1 ? 0.0 : -sign(r)

#==== M3 ====#

@ddffunc M3

M31(r) = 0.75-r^2
M32(r) = 0.5*(1.5-r)^2
@inline ddf_m3(r::Real) = r > 1.5 ? 0.0 : r <= 0.5 ? M31(r) : M32(r)

dM31(r) = -2*r
dM32(r) = r-1.5
@inline ddf_dm3(r::Real) = abs(r) > 1.5 ? 0.0 : abs(r) <= 0.5 ? dM31(r) : sign(r)*dM32(abs(r))

#==== YANG3 ====#

@ddffunc Yang3

Yang31(r) = 17/48+sqrt(3)π/108+r/4-r^2/4+(1-2*r)/16*sqrt(-12*r^2+12*r+1)-sqrt(3)/12*asin(sqrt(3)/2*(2*r-1))
Yang32(r) = 55/48-sqrt(3)π/108-13*r/12+r^2/4+(2*r-3)/48*sqrt(-12*r^2+36*r-23)+sqrt(3)/36*asin(sqrt(3)/2*(2*r-3))

@inline ddf_yang3(r::Real) = r > 2.0 ? 0.0 : r <= 1.0 ? Yang31(r) : Yang32(r)

function Base.show(io::IO, ddf::DDF{C,OVERDX}) where {C<:DDFType,OVERDX}
    print(io, "Discrete delta function operator of type $C, with spacing $(1.0/OVERDX)")
end
