abstract type DDFType end

struct DDF{C <: DDFType,OVERDX} end

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
Discrete delta function operator of type Whirl.Fields.Roma, with spacing 1.0

julia> ddf(1)
0.16666666666666666

julia> ddf(-1)
0.16666666666666666

julia> ddf.([-1,0,1])
3-element Array{Float64,1}:
 0.166667
 0.666667
 0.166667
```
"""
function DDF(;dx::Real=1.0,ddftype=Roma)
  DDF{ddftype,1.0/dx}()
end

# This macro creates an abstract type in the DDFType system for the named ddftype
# It also creates functions to evaluate the ddf in 1-d,2-d or 3-d.
macro ddffunc(ddftype)
    fname = Symbol("ddf_",lowercase(string(ddftype)))
    return esc(quote
            abstract type $ddftype <: DDFType end
            (::DDF{$ddftype,OVERDX})(x::T) where {OVERDX,T <: Real} =
                  $fname(abs(x)*OVERDX)*OVERDX
            (::DDF{$ddftype,OVERDX})(x::T,y::T) where {OVERDX,T <: Real} =
                          $fname(abs(x)*OVERDX)*OVERDX*$fname(abs(y)*OVERDX)*OVERDX
            (::DDF{$ddftype,OVERDX})(x::T,y::T,z::T) where {OVERDX,T <: Real} =
                          $fname(abs(x)*OVERDX)*OVERDX*$fname(abs(y)*OVERDX)*OVERDX*$fname(abs(z)*OVERDX)*OVERDX
            end)
end

@ddffunc Roma

roma1(r) = (1+sqrt(-3r^2+1))/3
roma2(r) = (5-3r-sqrt(1-3(1-r)^2))/6

ddf_roma(r::Real) = r > 1.5 ? 0.0 : r <= 0.5 ? roma1(r) : roma2(r)


@ddffunc Goza

ddf_goza(r::Real) = r > 14 ? 0.0 : exp(-π^2/36 * r^2)*sqrt(π/36)


@ddffunc Witchhat

ddf_witchhat(r::Real) = r > 1 ? 0.0 : 1-r


function Base.show(io::IO, ddf::DDF{C,OVERDX}) where {C<:DDFType,OVERDX}
    print(io, "Discrete delta function operator of type $C, with spacing $(1.0/OVERDX)")
end
