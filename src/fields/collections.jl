# Collections of data

# Collections of ScalarGridData

"""
    Edges

`Edges` is a wrapper for vector-valued data that lie at the faces of either dual cells or
primary cells. `Edges` type data have fields `u` and `v` for the components of the
vector field. These are the normal components of the vector field on the vertical
and horizontal faces of the corresponding cell.

# Constructors
- `Edges(C,dims)` creates a vector field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid.
- `Edges(C,w)` performs the same construction, but uses existing field data `w`
  of `GridData` type to determine the size of the grid.
-  Adding the `dtype=` keyword allows the data type of the field data to be
  changed. The default is `Float64`, but can be changed to, e.g., `ComplexF64`
"""
struct Edges{C <: CellType, NX, NY, T <: Number, DT} <: VectorGridData{NX,NY,T}
    data::DT
    u::XEdges{C,NX,NY,T}
    v::YEdges{C,NX,NY,T}
end

# Based on number of dual nodes, return the number of edges
edge_inds(::Type{C},   dualnodedims) where {C <: CellType} =
            xedge_inds(C,dualnodedims), yedge_inds(C,dualnodedims)


function (::Type{Edges{C,NX,NY,T,DT}})(data::AbstractVector{R}) where {C<: CellType,NX,NY,T<:Number,DT,R}
    udims = xedge_inds(C, (NX, NY))
    vdims = yedge_inds(C, (NX, NY))
    n0 = 0
    dims = udims
    dn = prod(dims)
    u = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = vdims
    dn = prod(dims)
    v = reshape(view(data,n0+1:n0+dn),dims)
    Edges{C, NX, NY,R,typeof(data)}(data, XEdges{C,NX,NY,R,typeof(u)}(u),
                                          YEdges{C,NX,NY,R,typeof(v)}(v))
end

function Edges(::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    udims = xedge_inds(C, dualnodedims)
    vdims = yedge_inds(C, dualnodedims)
    data = zeros(dtype,prod(udims)+prod(vdims))
    Edges{C,dualnodedims...,dtype,typeof(data)}(data)
end

@griddata(Edges,1)

#Base.IndexStyle(::Type{<:Edges}) = IndexLinear() # necessary?

function Base.show(io::IO, edges::Edges{C, NX, NY, T, DT}) where {C, NX, NY, T, DT}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(edges.u,1)), ny = $(size(edges.u,2)))"
    vdims = "(nx = $(size(edges.v,1)), ny = $(size(edges.v,2)))"
    println(io, "$C edges for a $nodedims cell grid of type $T data")
    println(io, "  Internal u-faces: $udims")
    print(io, "  Internal v-faces: $vdims")
end

function Base.show(io::IO, m::MIME"text/plain", edges::Edges)
    println(io,"$(typeof(edges)) data")
    println(io,"u (in grid orientation)")
    show(io,m,reverse(transpose(edges.u),dims=1))
    println(io)
    println(io,"v (in grid orientation)")
    show(io,m,reverse(transpose(edges.v),dims=1))
end

### Edge gradient ###

"""
    EdgeGradient

`EdgeGradient` is a wrapper for tensor-valued data that lie partly at the nodes of dual cells and
primary cells. `EdgeGradient` type data have fields `dudx`, `dudy`, `dvdx`, `dvdy` for the components of the
tensor field. The diagonal components lie at one set of nodes (e.g. Primal),
and the offdiagonal at the other set (e.g. Dual).

# Constructors
- `EdgeGradient(C,dims)` creates a tensor field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid.
- `EdgeGradient(C,w)` performs the same construction, but uses existing field data `w`
  of `GridData` type to determine the size of the grid.
-  Adding the `dtype=` keyword allows the data type of the field data to be
  changed. The default is `Float64`, but can be changed to, e.g., `ComplexF64`
"""
struct EdgeGradient{C <: CellType,D <: CellType, NX,NY, T<: Number, DT} <: GridData{NX,NY,T}
  data :: DT
  dudx :: Nodes{C,NX,NY,T}
  dvdy :: Nodes{C,NX,NY,T}
  dudy :: Nodes{D,NX,NY,T}
  dvdx :: Nodes{D,NX,NY,T}
end


# Should be able to clean these up...

function (::Type{EdgeGradient{C,D,NX,NY,T,DT}})(data::AbstractVector{R}) where {C<: CellType,D<:CellType,NX,NY,T<:Number,DT,R}
    dudxdims = dvdydims = node_inds(C, (NX,NY))
    dudydims = dvdxdims = node_inds(D, (NX,NY))

    n0 = 0
    dims = dudxdims
    dn = prod(dims)
    dudx = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = dudydims
    dn = prod(dims)
    dudy = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = dvdxdims
    dn = prod(dims)
    dvdx = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = dvdydims
    dn = prod(dims)
    dvdy = reshape(view(data,n0+1:n0+dn),dims)
    EdgeGradient{C,D,NX, NY,R,typeof(data)}(data, Nodes{C,NX,NY,R,typeof(dudx)}(dudx),
                                                  Nodes{C,NX,NY,R,typeof(dvdy)}(dvdy),
                                                  Nodes{D,NX,NY,R,typeof(dudy)}(dudy),
                                                  Nodes{D,NX,NY,R,typeof(dvdx)}(dvdx))
end

function EdgeGradient(::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    dudxdims = dvdydims = node_inds(C, dualnodedims)
    dudydims = dvdxdims = node_inds(othertype(C), dualnodedims)
    data = zeros(dtype,prod(dudxdims)+prod(dudydims)+prod(dvdxdims)+prod(dvdydims))
    EdgeGradient{C,othertype(C),dualnodedims...,dtype,typeof(data)}(data)
end

@griddata(EdgeGradient,2)

#Base.IndexStyle(::Type{<:EdgeGradient}) = IndexLinear()

function Base.show(io::IO, nodes::EdgeGradient{R, S, NX, NY, T, DT}) where {R, S, NX, NY, T, DT}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.dudx,1)), ny = $(size(nodes.dudx,2)))"
    vdims = "(nx = $(size(nodes.dudy,1)), ny = $(size(nodes.dudy,2)))"
    println(io, "Edge gradient data in a $nodedims cell grid of type $T data")
    println(io, "  du/dx and dv/dy as $R nodes : $udims")
    println(io, "  du/dy and dv/dx as $S nodes : $vdims")
end

function Base.show(io::IO, m::MIME"text/plain", nodes::EdgeGradient)
    println(io,"$(typeof(nodes)) data")
    println(io,"du/dx (in grid orientation)")
    show(io,m,reverse(transpose(nodes.dudx),dims=1))
    println(io)
    println(io,"du/dy (in grid orientation)")
    show(io,m,reverse(transpose(nodes.dudy),dims=1))
    println(io)
    println(io,"dv/dx (in grid orientation)")
    show(io,m,reverse(transpose(nodes.dvdx),dims=1))
    println(io)
    println(io,"dv/dy (in grid orientation)")
    show(io,m,reverse(transpose(nodes.dvdy),dims=1))
end

### Node pair ###

"""
    NodePair

`NodePair` is a wrapper for vector-valued data that lie at the nodes of dual cells and
primal cells. `NodePair` type data have fields `u` and `v` for the components of the
vector field. These are the normal components of a vector field on nodes that
form the faces of a virtual cell centered at one of the faces of the primal cell.

# Constructors
- `NodePair(C,dims)` creates a vector field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid.
- `NodePair(C,w)` performs the same construction, but uses existing field data `w`
  of `GridData` type to determine the size of the grid.
-  Adding the `dtype=` keyword allows the data type of the field data to be
  changed. The default is `Float64`, but can be changed to, e.g., `ComplexF64`
"""
struct NodePair{C <: CellType,D <: CellType, NX,NY,T, DT} <: GridData{NX,NY,T}
  data :: DT
  u :: Nodes{C,NX,NY,T}
  v :: Nodes{D,NX,NY,T}
end

function (::Type{NodePair{C,D,NX,NY,T,DT}})(data::AbstractVector{R}) where {C<: CellType,D<: CellType,NX,NY,T<:Number,DT,R}
    udims = node_inds(C, (NX,NY))
    vdims = node_inds(D, (NX,NY))
    n0 = 0
    dims = udims
    dn = prod(dims)
    u = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = vdims
    dn = prod(dims)
    v = reshape(view(data,n0+1:n0+dn),dims)
    NodePair{C, D, NX, NY,R,typeof(data)}(data, Nodes{C,NX,NY,R,typeof(u)}(u),
                                                Nodes{D,NX,NY,R,typeof(v)}(v))
end

function NodePair(::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    udims = node_inds(C, dualnodedims)
    vdims = node_inds(othertype(C), dualnodedims)
    data = zeros(dtype,prod(udims)+prod(vdims))
    NodePair{C,othertype(C),dualnodedims...,dtype,typeof(data)}(data)
end

@griddata(NodePair,2)

#Base.IndexStyle(::Type{<:NodePair}) = IndexLinear() # necessary?

function Base.show(io::IO, nodes::NodePair{R, S, NX, NY, T}) where {R, S, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.u,1)), ny = $(size(nodes.u,2)))"
    vdims = "(nx = $(size(nodes.v,1)), ny = $(size(nodes.v,2)))"
    println(io, "($R,$S) node data pair in a $nodedims cell grid of type $T data")
    println(io, "  Number of $R nodes: $udims")
    println(io, "  Number of $S nodes: $vdims")
end

function Base.show(io::IO, m::MIME"text/plain", nodes::NodePair)
    println(io,"$(typeof(nodes)) data")
    println(io,"u (in grid orientation)")
    show(io,m,reverse(transpose(nodes.u),dims=1))
    println(io)
    println(io,"v (in grid orientation)")
    show(io,m,reverse(transpose(nodes.v),dims=1))
end
