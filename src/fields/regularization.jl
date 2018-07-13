
struct Regularize{N,DV,F}

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

  "Symmetry flag"
  _issymmetric :: Bool

end


"""
    Regularize(x,y,dx,[ddftype=Roma],[I0=(1,1)], [weights=1.0], [filter=false],
                       [issymmetric=false])

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

If the optional Boolean parameter `issymmetric` is set to `true`, then the
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
                    issymmetric::Bool = false) where {T<:Real}
  n = length(x)
  @assert length(y)==n
  if !issymmetric
    if typeof(weights) == T
      wtvec = similar(x)
      fill!(wtvec,weights)
    else
      @assert length(weights)==n
      wtvec = deepcopy(weights)
    end
  else
    # if the regularization and interpolation are symmetric, then the
    # weights are automatically set to be the cell area in order to cancel it
    # in the denominator of the regularization operator.
    wtvec = similar(x)
    fill!(wtvec,dx*dx)
  end

  Regularize{length(x),dx*dx,filter}(x/dx+I0[1],y/dx+I0[2],
                      wtvec,zeros(T,n),DDF(ddftype=ddftype,dx=1.0),issymmetric)
end

Regularize(x::T,y::T,a...;b...) where {T<:Real} = Regularize([x],[y],a...;b...)

Regularize(x::VectorData,a...;b...) = Regularize(x.u,x.v,a...;b...)

function Base.show(io::IO, H::Regularize{N,DV,F}) where {N,DV,F}
    filter = F ? "filtered" : "non-filtered"
    op = H._issymmetric ? "Symmetric regularization/interpolation" : "Regularization/interpolation"
    println(io, "$op operator with $filter interpolation")
    println(io, "  $N points in grid with cell area $(sprint(showcompact,DV))")
end

"""
    RegularizationMatrix(H::Regularize,f::ScalarData/VectorData,u::Nodes/Edges) -> Hmat

Construct and store a matrix representation of regularization associated with `H`
for data of type `f` to data of type `u`. The resulting matrix `Hmat` can then be
used to apply on point data of type `f` to regularize it to grid data of type `u`,
using `A_mul_B!(u,Hmat,f)`. It can also be used as just `Hmat*f`.

If `H` is a symmetric regularization and interpolation operator, then this
actually returns a tuple `Hmat, Emat`, where `Emat` is the interpolation matrix.
"""
struct RegularizationMatrix{TU,TF} <: AbstractMatrix{Float64}
  M :: SparseMatrixCSC{Float64,Int64}
end


"""
    InterpolationMatrix(H::Regularize,u::Nodes/Edges,f::ScalarData/VectorData) -> Emat

Construct and store a matrix representation of interpolation associated with `H`
for data of type `u` to data of type `f`. The resulting matrix `Emat` can then be
used to apply on grid data of type `u` to interpolate it to point data of type `f`,
using `A_mul_B!(f,Emat,u)`. It can also be used as just `Emat*u`.
"""
struct InterpolationMatrix{TU,TF} <: AbstractMatrix{Float64}
  M :: SparseMatrixCSC{Float64,Int64}
end

@wraparray RegularizationMatrix M
@wraparray InterpolationMatrix M


# Regularization and interpolation operators of vector data to edges
ftype = :(VectorData{N})
for (ctype,dunx,duny,dvnx,dvny,shiftux,shiftuy,shiftvx,shiftvy) in vectorlist

# Regularization
  @eval function (H::Regularize{N,DV,F})(target::$ctype,source::$ftype) where {N,DV,F,NX,NY}
        fill!(target.u,0.0)
        @inbounds for y in 1:NY-$duny, x in 1:NX-$dunx
          H.buffer .= H.ddf.(x-$shiftux-H.x,y-$shiftuy-H.y)
          target.u[x,y] = dot(H.buffer,source.u.*H.wgt)/DV
        end
        fill!(target.v,0.0)
        @inbounds for y in 1:NY-$dvny, x in 1:NX-$dvnx
          H.buffer .= H.ddf.(x-$shiftvx-H.x,y-$shiftvy-H.y)
          target.v[x,y] = dot(H.buffer,source.v.*H.wgt)/DV
        end
        target
  end

# Interpolation
  @eval function (H::Regularize{N,DV,false})(target::$ftype,
                                       source::$ctype) where {N,DV,NX,NY}
    target.u .= target.v .= zeros(Float64,N)
    @inbounds for y in 1:NY-$duny, x in 1:NX-$dunx
      H.buffer .= H.ddf.(x-$shiftux-H.x,y-$shiftuy-H.y)
      target.u .+= H.buffer*source.u[x,y]
    end
    @inbounds for y in 1:NY-$dvny, x in 1:NX-$dvnx
      H.buffer .= H.ddf.(x-$shiftvx-H.x,y-$shiftvy-H.y)
      target.v .+= H.buffer*source.v[x,y]
    end
    target
  end

# Interpolation with filtering
  @eval function (H::Regularize{N,DV,true})(target::$ftype,
                                      source::$ctype) where {N,DV,NX,NY}
    target.u .= target.v .= zeros(Float64,N)
    @inbounds for y in 1:NY-$duny, x in 1:NX-$dunx
      H.buffer .= H.ddf.(x-$shiftux-H.x,y-$shiftuy-H.y)
      w = dot(H.buffer,H.wgt)/DV
      w = w ≢ 0.0 ? source.u[x,y]/w : 0.0
      target.u .+= H.buffer*w
    end
    @inbounds for y in 1:NY-$dvny, x in 1:NX-$dvnx
      H.buffer .= H.ddf.(x-$shiftvx-H.x,y-$shiftvy-H.y)
      w = dot(H.buffer,H.wgt)/DV
      w = w ≢ 0.0 ? source.v[x,y]/w : 0.0
      target.v .+= H.buffer*w
    end
    target
  end

  # Construct regularization matrix
  @eval function RegularizationMatrix(H::Regularize{N,DV,F},src::$ftype,target::$ctype) where {N,DV,F,NX,NY}

    #Hmat = (spzeros(length(target.u),length(src.u)),spzeros(length(target.v),length(src.v)))
    lenu = length(target.u)
    lenv = length(target.v)
    Hmat = spzeros(lenu+lenv,2N)
    g = deepcopy(src)
    v = deepcopy(target)
    g.u .= g.v .= zeros(Float64,N)
    for i = 1:N
      g.u[i] = 1.0
      g.v[i] = 1.0
      H(v,g)
      Hmat[1:lenu,i]           = sparsevec(v.u)
      Hmat[lenu+1:lenu+lenv,i+N] = sparsevec(v.v)
      g.u[i] = 0.0
      g.v[i] = 0.0
    end
    if H._issymmetric
      # In symmetric case, these matrices are identical. (Interpolation is stored
      # as its transpose.)
      return RegularizationMatrix{$ctype,$ftype}(Hmat),InterpolationMatrix{$ctype,$ftype}(Hmat)
    else
      return RegularizationMatrix{$ctype,$ftype}(Hmat)
    end
  end

  # Construct interpolation matrix
  @eval function InterpolationMatrix(H::Regularize{N,DV,false},src::$ctype,target::$ftype) where {N,DV,NX,NY}

    # note that we store interpolation matrices in the same shape as regularization matrices
    #Emat = (spzeros(length(src.u),length(target.u)),spzeros(length(src.v),length(target.v)))
    lenu = length(src.u)
    lenv = length(src.v)
    Emat = spzeros(lenu+lenv,2N)
    g = deepcopy(target)
    v = deepcopy(src)
    g.u .= g.v .= zeros(Float64,N)
    for i = 1:N
      g.u[i] = DV/H.wgt[i]  # unscale for interpolation
      g.v[i] = DV/H.wgt[i]  # unscale for interpolation
      H(v,g)
      Emat[1:lenu,i]           = sparsevec(v.u)
      Emat[lenu+1:lenu+lenv,i+N] = sparsevec(v.v)
      g.u[i] = 0.0
      g.v[i] = 0.0
    end
    InterpolationMatrix{$ctype,$ftype}(Emat)
  end

  # Construct interpolation matrix with filtering
  @eval function InterpolationMatrix(H::Regularize{N,DV,true},src::$ctype,target::$ftype) where {N,DV,NX,NY}

    # note that we store interpolation matrices in the same shape as regularization matrices
    #Emat = (spzeros(length(src.u),length(target.u)),spzeros(length(src.v),length(target.v)))
    lenu = length(src.u)
    lenv = length(src.v)
    Emat = spzeros(lenu+lenv,2N)
    g = deepcopy(target)
    v = deepcopy(src)
    fill!(g,1.0)
    H(v,g)
    wtu = sparsevec(v.u)
    wtu.nzval .= 1./wtu.nzval
    wtv = sparsevec(v.v)
    wtv.nzval .= 1./wtv.nzval
    fill!(g,0.0)
    for i = 1:N
      g.u[i] = DV/H.wgt[i]  # unscale for interpolation
      g.v[i] = DV/H.wgt[i]  # unscale for interpolation
      H(v,g)
      Emat[1:lenu,i]             = wtu.*sparsevec(v.u)
      Emat[lenu+1:lenu+lenv,i+N] = wtv.*sparsevec(v.v)
      g.u[i] = 0.0
      g.v[i] = 0.0
    end
    InterpolationMatrix{$ctype,$ftype}(Emat)
  end

  @eval function A_mul_B!(u::$ctype,Hmat::RegularizationMatrix{$ctype,$ftype},f::$ftype) where {NX,NY,N}
    fill!(u,0.0)
    I,J,V = findnz(Hmat.M)
    for (cnt,v) in enumerate(V)
      u[I[cnt]] += v*f[J[cnt]]
    end
    u
  end

  @eval function A_mul_B!(f::$ftype,Emat::InterpolationMatrix{$ctype,$ftype},u::$ctype) where {NX,NY,N}
    fill!(f,0.0)
    I,J,V = findnz(Emat.M)
    for (cnt,v) in enumerate(V)
      f[J[cnt]] .+= v*u[I[cnt]]
    end
    f
  end

end


# Nodal type
ftype = :(ScalarData{N})
for (ctype,dnx,dny,shiftx,shifty) in scalarlist

# Regularization
  @eval function (H::Regularize{N,DV,F})(target::$ctype,source::$ftype) where {N,DV,F,NX,NY}
    fill!(target,0.0)
    @inbounds for y in 1:NY-$dny, x in 1:NX-$dnx
      H.buffer .= H.ddf.(x-$shiftx-H.x,y-$shifty-H.y)
      target[x,y] = dot(H.buffer,source.data.*H.wgt)/DV
    end
    target
  end

# Interpolation
  @eval function (H::Regularize{N,DV,false})(target::$ftype,
                                       source::$ctype) where {N,DV,NX,NY}
    target .= zeros(Float64,N)
    @inbounds for y in 1:NY-$dny, x in 1:NX-$dnx
      H.buffer .= H.ddf.(x-$shiftx-H.x,y-$shifty-H.y)
      target .+= H.buffer*source[x,y]
    end
    target
  end

# Interpolation with filtering
  @eval function (H::Regularize{N,DV,true})(target::$ftype,
                                       source::$ctype) where {N,DV,NX,NY}
    target .= zeros(Float64,N)
    @inbounds for y in 1:NY-$dny, x in 1:NX-$dnx
      H.buffer .= H.ddf.(x-$shiftx-H.x,y-$shifty-H.y)
      w = dot(H.buffer,H.wgt)/DV
      w = w ≢ 0.0 ? source[x,y]/w : 0.0
      target .+= H.buffer*w
    end
    target
  end

  # Construct regularization matrix
  @eval function RegularizationMatrix(H::Regularize{N,DV,F},
    f::$ftype,
    u::$ctype) where {N,DV,F,NX,NY}

    Hmat = spzeros(length(u),length(f))
    g = deepcopy(f)
    v = deepcopy(u)
    fill!(g,0.0)
    for i = 1:N
      g[i] = 1.0
      Hmat[:,i] = sparsevec(H(v,g))
      g[i] = 0.0
    end
    if H._issymmetric
      # In symmetric case, these matrices are identical. (Interpolation is stored
      # as its transpose.)
      return RegularizationMatrix{$ctype,$ftype}(Hmat),InterpolationMatrix{$ctype,$ftype}(Hmat)
    else
      return RegularizationMatrix{$ctype,$ftype}(Hmat)
    end
  end

  # Construct interpolation matrix
  @eval function InterpolationMatrix(H::Regularize{N,DV,false},
    u::$ctype,
    f::$ftype) where {N,DV,NX,NY}

    Emat = spzeros(length(u),length(f))
    g = deepcopy(f)
    v = deepcopy(u)
    fill!(g,0.0)
    for i = 1:N
      g[i] = DV/H.wgt[i]  # unscale for interpolation
      Emat[:,i] = sparsevec(H(v,g))
      g[i] = 0.0
    end
    InterpolationMatrix{$ctype,$ftype}(Emat)
  end

  # Construct interpolation matrix with filtering
  @eval function InterpolationMatrix(H::Regularize{N,DV,true},
    u::$ctype,
    f::$ftype) where {N,DV,NX,NY}

    Emat = spzeros(length(u),length(f))
    g = deepcopy(f)
    v = deepcopy(u)
    fill!(g,1.0)
    wt = sparsevec(H(v,g))
    wt.nzval .= 1./wt.nzval
    fill!(g,0.0)
    for i = 1:N
      g[i] = DV/H.wgt[i]  # unscale for interpolation
      Emat[:,i] = wt.*sparsevec(H(v,g))
      g[i] = 0.0
    end
    InterpolationMatrix{$ctype,$ftype}(Emat)
  end

  @eval function A_mul_B!(u::$ctype,Hmat::RegularizationMatrix{$ctype,$ftype},f::$ftype) where {NX,NY,N}
    fill!(u,0.0)
    I,J,V = findnz(Hmat.M)
    for (cnt,v) in enumerate(V)
      u[I[cnt]] += v*f[J[cnt]]
    end
    u

  end

  @eval function A_mul_B!(f::$ftype,Emat::InterpolationMatrix{$ctype,$ftype},u::$ctype) where {NX,NY,N}
    fill!(f,0.0)
    I,J,V = findnz(Emat.M)
    for (cnt,v) in enumerate(V)
      f[J[cnt]] += v*u[I[cnt]]
    end
    f

  end


end

(*)(Hmat::RegularizationMatrix{TU,TF},src::TF) where {TU,TF<:Union{ScalarData,VectorData}} =
        A_mul_B!(TU(),Hmat,src)

(*)(Emat::InterpolationMatrix{TU,TF},src::TU) where {TU<:Union{Nodes,Edges},TF<:Union{ScalarData,VectorData}} =
                A_mul_B!(TF(),Emat,src)


function Base.show(io::IO, H::RegularizationMatrix{TU,TF}) where {TU,TF}
    print(io, "Regularization matrix acting on type $TF and returning type $TU")
end

function Base.show(io::IO, H::InterpolationMatrix{TU,TF}) where {TU,TF}
    print(io, "Interpolation matrix acting on type $TU and returning type $TF")
end
