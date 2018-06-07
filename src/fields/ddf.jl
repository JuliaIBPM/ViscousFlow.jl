abstract type DDFType end

struct DDF{C <: DDFType,DX} end

function DDF(;dx::Real=1.0,ddftype=Roma)
  DDF{ddftype,dx}()
end

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

@ddffunc Tophat

tophat(r) = 1-r;

function ddf_tophat(r::Real)
  rr = abs(r)
  return rr > 1 ? 0.0 : tophat(rr)
end
