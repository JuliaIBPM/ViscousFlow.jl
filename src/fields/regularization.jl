
struct Regularize{N,DV,S}

  "x values of points, normalized to grid index space"
  x :: Vector{Float64}

  "y values of points, normalized to grid index space"
  y :: Vector{Float64}

  "weights for each point (e.g. arclengths)"
  wgt :: Vector{Float64}

  "buffer space"
  buffer :: Vector{Float64}

  "Discrete Delta function"
  ddf :: DDF
end


"""
    Regularize(x,y,dx,[ddftype=Roma],[I0=(1,1)], [weights=1.0], [filter=false],
                       [symmetric=false])

Constructor to set up an operator for regularizing and interpolating data from/to
points immersed in the grid to/from fields on the grid itself. The supplied
`x` and `y` represent physical coordinates of the immersed points, and `dx`
denotes a uniform physical cell size of the grid. The separate arguments `x` and
`y` can be replaced by a single argument `X` of type `VectorData` holding the
coordinates.

The operations of regularization and interpolation are carried out with a discrete
delta function (ddf), which defaults to the type `Roma`. Others are also possible,
such as `Goza`. The optional tuple
`I0` represents the indices of the primary node that coincides with `(x,y) = (0,0)`.
This defaults to `(1,1)`, which leaves one layer of ghost (dual) cells and sets
the physical origin in the lower left corner of the grid of interior dual cells.

Another optional parameter, `weights`, sets the weight of each point in the
regularization. This would generally be set with, say, the differential arc
length for regularization of data on a curve. It can be a vector (of the same length
as x and y) or a scalar if uniform. It defaults to 1.0.

The optional Boolean parameter `filter` can be set to `true` if it is desired to
apply filtering (see Goza et al, J Comput Phys 2016) to the grid data before
interpolating. This is generally only used in the context of preconditioning
the solution for forces on the immersed points.

If the optional Boolean parameter `symmetric` is set to `true`, then the
regularization and interpolation are constructed to be transposes of each other.
Note that this option overrides any supplied weights. The default of this
parameter is `false`.

The resulting operator can be used in either direction, regularization and
interpolation, with the first argument representing the *target* (the entity
to regularize/interpolate to), and the second argument
the *source* (the entity to regularize/interpolate from). The regularization
does not use the filtering option.

# Example

In the example below, we set up a 12 x 12 grid. Using the default value for `I0`
and setting `dx = 0.1`, the physical dimensions of the non-ghost part of the grid
are 1.0 x 1.0. Three points are set up in the interior, and a vector field is assigned
to them, with the x component of each of them set to 1.0. These data are regularized
to a field of primal edges on the grid.

```jldoctest
julia> x = [0.25,0.75,0.25]; y = [0.75,0.25,0.25];

julia> X = VectorData(x,y);

julia> q = Edges(Primal,(12,12));

julia> dx = 0.1;

julia> H = Regularize(x,y,dx)
Regularization/interpolation operator with non-filtered interpolation
  3 points in grid with cell area 0.01

julia> f = VectorData(X);

julia> fill!(f.u,1.0);

julia> H(q,f)
Whirl.Fields.Edges{Whirl.Fields.Primal,12,12} data
u (in grid orientation):
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
v (in grid orientation):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function Regularize(x::Vector{T},y::Vector{T},dx::T;
                    ddftype=Roma,
                    I0::Tuple{Int,Int}=(1,1),
                    weights::Union{T,Vector{T}}=1.0,
                    filter::Bool = false,
                    symmetric::Bool = false) where {T<:Real}
  n = length(x)
  @assert length(y)==n
  if !symmetric
    if typeof(weights) == T
      wtvec = similar(x)
      fill!(wtvec,weights)
    else
      @assert length(weights)==n
      wtvec = deepcopy(weights)
    end
  else
    wtvec = similar(x)
    fill!(wtvec,dx*dx)
  end

  Regularize{length(x),dx*dx,filter}(x/dx+I0[1],y/dx+I0[2],wtvec,zeros(T,n),DDF(ddftype=ddftype,dx=1.0))
end

Regularize(x::T,y::T,a...;b...) where {T<:Real} = Regularize([x],[y],a...;b...)

Regularize(x::VectorData,a...;b...) = Regularize(x.u,x.v,a...;b...)

function Base.show(io::IO, H::Regularize{N,DV,S}) where {N,DV,S}
    filter = S ? "filtered" : "non-filtered"
    println(io, "Regularization/interpolation operator with $filter interpolation")
    println(io, "  $N points in grid with cell area $(sprint(showcompact,DV))")
end

# These regularization operations should be easy to macro-generate

function (H::Regularize{N,DV})(target::Edges{Primal,NX,NY},source::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    target.u[x,y] = dot(H.buffer,source.u.*H.wgt)/DV
  end
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    target.v[x,y] = dot(H.buffer,source.v.*H.wgt)/DV
  end
  target
end

function (H::Regularize{N,DV})(target::Edges{Dual,NX,NY},source::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    target.u[x,y] = dot(H.buffer,source.u.*H.wgt)/DV
  end
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    target.v[x,y] = dot(H.buffer,source.v.*H.wgt)/DV
  end
  target
end

function (H::Regularize{N,DV})(target::Tuple{Nodes{Dual,NX,NY},Nodes{Primal,NX,NY}},source::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target[1][x,y] = dot(H.buffer,source.u.*H.wgt)/DV
  end
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target[2][x,y] = dot(H.buffer,source.v.*H.wgt)/DV
  end
  target
end

function (H::Regularize{N,DV})(target::Tuple{Nodes{Primal,NX,NY},Nodes{Dual,NX,NY}},source::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target[1][x,y] = dot(H.buffer,source.u.*H.wgt)/DV
  end
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target[2][x,y] = dot(H.buffer,source.v.*H.wgt)/DV
  end
  target
end

function (H::Regularize{N,DV})(target::Nodes{Primal,NX,NY},source::ScalarData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target[x,y] = dot(H.buffer,source.data.*H.wgt)/DV
  end
  target
end

function (H::Regularize{N,DV})(target::Nodes{Dual,NX,NY},source::ScalarData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target[x,y] = dot(H.buffer,source.data.*H.wgt)/DV
  end
  target
end

# Interpolation functions

function (H::Regularize{N,DV,false})(target::VectorData{N},
                                     source::Edges{Primal,NX,NY}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    target.u .+= H.buffer*source.u[x,y]
  end
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    target.v .+= H.buffer*source.v[x,y]
  end
  target
end

function (H::Regularize{N,DV,false})(target::VectorData{N},
                                     source::Edges{Dual,NX,NY}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    target.u .+= H.buffer*source.u[x,y]
  end
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    target.v .+= H.buffer*source.v[x,y]
  end
  target
end

function (H::Regularize{N,DV,false})(target::VectorData{N},
                                     source::Tuple{Nodes{Dual,NX,NY},Nodes{Primal,NX,NY}}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target.u .+= H.buffer*source[1][x,y]
  end
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target.v .+= H.buffer*source[2][x,y]
  end
  target
end

function (H::Regularize{N,DV,false})(target::VectorData{N},
                                     source::Tuple{Nodes{Primal,NX,NY},Nodes{Dual,NX,NY}}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target.u .+= H.buffer*source[1][x,y]
  end
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target.v .+= H.buffer*source[2][x,y]
  end
  target
end


function (H::Regularize{N,DV,false})(target::ScalarData{N},
                                     source::Nodes{Primal,NX,NY}) where {N,DV,NX,NY}
  target .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    target .+= H.buffer*source[x,y]
  end
  target
end

function (H::Regularize{N,DV,false})(target::ScalarData{N},
                                     source::Nodes{Dual,NX,NY}) where {N,DV,NX,NY}
  target .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    target .+= H.buffer*source[x,y]
  end
  target
end

# Interpolation with filtering

function (H::Regularize{N,DV,true})(target::VectorData{N},
                                    source::Edges{Primal,NX,NY}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source.u[x,y]/w : 0.0
    target.u .+= H.buffer*w
  end
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source.v[x,y]/w : 0.0
    target.v .+= H.buffer*w
  end
  target
end

function (H::Regularize{N,DV,true})(target::VectorData{N},
                                    source::Edges{Dual,NX,NY}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-0.5-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source.u[x,y]/w : 0.0
    target.u .+= H.buffer*w
  end
  @inbounds for y in 1:NY-1, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source.v[x,y]/w : 0.0
    target.v .+= H.buffer*w
  end
  target
end

function (H::Regularize{N,DV,true})(target::VectorData{N},
                                     source::Tuple{Nodes{Dual,NX,NY},Nodes{Primal,NX,NY}}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[1][x,y]/w : 0.0
    target.u .+= H.buffer*w
  end
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[2][x,y]/w : 0.0
    target.v .+= H.buffer*w
  end
  target
end

function (H::Regularize{N,DV,true})(target::VectorData{N},
                                     source::Tuple{Nodes{Primal,NX,NY},Nodes{Dual,NX,NY}}) where {N,DV,NX,NY}
  target.u .= target.v .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[1][x,y]/w : 0.0
    target.u .+= H.buffer*w
  end
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[2][x,y]/w : 0.0
    target.v .+= H.buffer*w
  end
  target
end


function (H::Regularize{N,DV,true})(target::ScalarData{N},
                                     source::Nodes{Primal,NX,NY}) where {N,DV,NX,NY}
  target .= zeros(Float64,N)
  @inbounds for y in 1:NY-1, x in 1:NX-1
    H.buffer .= H.ddf.(x-H.x,y-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[x,y]/w : 0.0
    target .+= H.buffer*w
  end
  target
end

function (H::Regularize{N,DV,true})(target::ScalarData{N},
                                    source::Nodes{Dual,NX,NY}) where {N,DV,NX,NY}
  target .= zeros(Float64,N)
  @inbounds for y in 1:NY, x in 1:NX
    H.buffer .= H.ddf.(x-0.5-H.x,y-0.5-H.y)
    w = dot(H.buffer,H.wgt)/DV
    w = w ≢ 0.0 ? source[x,y]/w : 0.0
    target .+= H.buffer*w
  end
  target
end
