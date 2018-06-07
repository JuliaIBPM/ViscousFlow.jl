abstract type DDFType end

struct DDF{C <: DDFType,DX} end

"""
    DDF([ddftype=Fields.Roma],[dx=1.0])

Construct a discrete delta function operator. This is generally only needed
internally by the `Regularize` operator, so the user doesn't have much need
for accessing this directly. The default DDF is the `Roma` function, which
has a support of 3 grid cells. Other choices are the `Goza` operator, which is
a truncated Gaussian with 28 cells support, and the `Witchhat`, which has 2
cells support. The resulting operator is evaluated with one, two or three
coordinates, producing, respectively, 1-d, 2-d, and 3-d smeared delta functions.
It can also be called with the usual Julia vectorized dot notation with arrays of
arguments.
The optional cell spacing argument `dx` rescales the coordinates by this spacing,
and the result is also rescaled by this spacing (raised to the number of dimensions).

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
  DDF{ddftype,dx}()
end

# This macro creates an abstract type in the DDFType system for the named ddftype
# It also creates functions to evaluate the ddf in 1-d,2-d or 3-d.
macro ddffunc(ddftype)
    fname = Symbol("ddf_",lowercase(string(ddftype)))
    return esc(quote
            abstract type $ddftype <: DDFType end
            (::DDF{$ddftype,DX})(x::T) where {DX,T <: Real} =
                  $fname(abs(x)/DX)/DX
            (::DDF{$ddftype,DX})(x::T,y::T) where {DX,T <: Real} =
                          $fname(abs(x)/DX)/DX*$fname(abs(y)/DX)/DX
            (::DDF{$ddftype,DX})(x::T,y::T,z::T) where {DX,T <: Real} =
                          $fname(abs(x)/DX)/DX*$fname(abs(y)/DX)/DX*$fname(abs(z)/DX)/DX
            #=
            (::DDF{$ddftype,DX})(x::AbstractVector{T}) where {DX,T <: Real} =
                          $fname.(abs.(x)/DX)
            (::DDF{$ddftype,DX})(x::AbstractVector{T},
                                 y::AbstractVector{T}) where {DX,T <: Real} =
                          $fname.(abs.(x)/DX).*$fname.(abs.(y)/DX)
            (::DDF{$ddftype,DX})(x::AbstractVector{T},
                                 y::AbstractVector{T},
                                 z::AbstractVector{T}) where {DX,T <: Real} =
                          $fname.(abs.(x)/DX).*$fname.(abs.(y)/DX).*$fname.(abs.(z)/DX)
            =#
            end)
end

@ddffunc Roma

roma1(r) = (1+sqrt(-3r^2+1))/3;
roma2(r) = (5-3r-sqrt(1-3(1-r)^2))/6;

function ddf_roma(r::Real)
    rr = abs(r)
    if rr > 1.5
      return 0.0
    end
    if rr <= 0.5
      return roma1(rr)
    else
      return roma2(rr)
    end
end



#=
function ddf_roma(r::Array{Float64,1})
  rr = abs.(r)

  f = zeros(rr)
  pts = find(x -> x<=0.5,rr)
  f[pts] = roma1.(rr[pts])

  pts = find(x -> 0.5 < x < 1.5,rr)
  f[pts] = roma2.(rr[pts])

  f

end
=#

@ddffunc Goza

goza(r) = exp(-π^2/36 * r^2)*sqrt(π/36);

function ddf_goza(r::Real)
  rr = abs(r)
  return rr > 14 ? 0.0 : goza(rr)
end

#=
function ddf_goza(r::Array{Float64,1})
  rr = abs.(r)

  f = zeros(rr)
  pts = find(x -> x<=14,rr)
  f[pts] = goza.(rr[pts])

  f
end
=#

@ddffunc Witchhat

witchhat(r) = 1-r;

function ddf_witchhat(r::Real)
  rr = abs(r)
  return rr > 1 ? 0.0 : witchhat(rr)
end

function Base.show(io::IO, ddf::DDF{C,DX}) where {C<:DDFType,DX}
    print(io, "Discrete delta function operator of type $C, with spacing $DX")
end
