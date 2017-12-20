module Process

import Whirl2d
import Whirl2d:@get
#import Whirl2d.Grids
#import Whirl2d.Bodies
#import Whirl2d.Systems

using Interpolations

struct SamplePoints
  x :: Vector{Float64}
  y :: Vector{Float64}
end

SamplePoint(x::Float64,y::Float64) = SamplePoints([x],[y])

function sample(pts::SamplePoints,u::Vector{Array{T,2}},x::Vector{Float64},y::Vector{Float64}) where {T}
  usamp = [T[] for x in pts.x]

  for (i,ui) in enumerate(u)
    uitp = interpolate((x,y),ui,Gridded(Linear()))
    for j = 1:length(usamp)
      push!(usamp[j],uitp[pts.x[j],pts.y[j]])
    end
  end
  return usamp

end

sample(pts::SamplePoints,u::Array{T,2},x::Vector{Float64},y::Vector{Float64}) where {T} =
      sample(pts,[u],x,y)



end
