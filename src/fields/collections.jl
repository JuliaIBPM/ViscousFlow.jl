# Collections of data

# To do
# - wrap a vector, as with Edges, and use view/reshape for the components
# - eliminate the second CellType parameter, since it is not really needed

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
    dudydims = dvdxdims = node_inds(othertype(C), (NX,NY))

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


# These will get handled by @griddata when the second type parameter is removed
(::Type{EdgeGradient{R,S,NX,NY,T,DT}})() where {R<:CellType,S<:CellType,NX,NY,T,DT} = EdgeGradient(R, (NX, NY),dtype=T)

(::Type{EdgeGradient{R,S,NX,NY,T}})() where {R<:CellType,S<:CellType,NX,NY,T} = EdgeGradient(R, (NX, NY),dtype=T)

EdgeGradient(C, ::GridData{NX,NY,T}; dtype=T) where {NX, NY,T} = EdgeGradient(C, (NX,NY),dtype=dtype)

Base.similar(::EdgeGradient{R,S,NX,NY,T,DT};element_type=T) where {R,S,NX,NY,T,DT} = EdgeGradient(R, (NX, NY),dtype=element_type)


# These should get handled same way as other GridData (with possible exception of ScalarGridData)
Base.size(A::EdgeGradient) = size(A.data)
@propagate_inbounds Base.getindex(A::EdgeGradient,i::Int) = getindex(A.data,i)
@propagate_inbounds Base.setindex!(A::EdgeGradient{R,S,NX,NY,T}, v, i::Int) where {R,S,NX,NY,T}= setindex!(A.data,convert(T,v),i)
Base.IndexStyle(::Type{<:EdgeGradient}) = IndexLinear()

Base.parent(A::EdgeGradient) = A.data
Base.parentindices(A::EdgeGradient) = parentindices(A.data)

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


struct NodePair{C <: CellType,D <: CellType, NX,NY,T, DT} <: GridData{NX,NY,T}
  data :: DT
  u :: Nodes{C,NX,NY,T}
  v :: Nodes{D,NX,NY,T}
end

function (::Type{NodePair{C,D,NX,NY,T,DT}})(data::AbstractVector{R}) where {C<: CellType,D<: CellType,NX,NY,T<:Number,DT,R}
    udims = node_inds(C, (NX,NY))
    vdims = node_inds(othertype(C), (NX,NY))
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

# These will get handled by @griddata when the second type parameter is removed
(::Type{NodePair{R,S,NX,NY,T,DT}})() where {R<:CellType,S<:CellType,NX,NY,T,DT} = NodePair(R, (NX, NY),dtype=T)

(::Type{NodePair{R,S,NX,NY,T}})() where {R<:CellType,S<:CellType,NX,NY,T} = NodePair(R, (NX, NY),dtype=T)

NodePair(C, ::GridData{NX,NY,T}; dtype=T) where {NX, NY,T} = NodePair(C, (NX,NY),dtype=dtype)

Base.similar(::NodePair{R,S,NX,NY,T,DT};element_type=T) where {R,S,NX,NY,T,DT} = NodePair(R, (NX, NY),dtype=element_type)

Base.size(A::NodePair) = size(A.data)
@propagate_inbounds Base.getindex(A::NodePair{C,D,NX,NY,T},i::Int) where {C,D,NX,NY,T} = getindex(A.data,i)
@propagate_inbounds Base.setindex!(A::NodePair{C,D,NX,NY,T}, v, i::Int) where {C,D,NX,NY,T} = setindex!(A.data,convert(T,v),i)
Base.IndexStyle(::Type{<:NodePair}) = IndexLinear() # necessary?


#=
function NodePair(R::Type{C},S::Type{D},dualnodedims::Tuple{Int,Int};dtype=Float64) where {C <: CellType, D <: CellType}
    udims = node_inds(R, dualnodedims)
    vdims = node_inds(S, dualnodedims)

    NodePair{R, S, dualnodedims...,dtype}(
            Nodes(R,dualnodedims,dtype=dtype),
            Nodes(S,dualnodedims,dtype=dtype)
            )
end
NodePair(R::Type{C}, U::Type{D},nodes::Nodes{S,NX,NY,T}) where
        {C <: CellType, D <: CellType, S <: CellType, NX,NY,T} = NodePair(R, U, (NX, NY),dtype=T)
NodePair(R::Type{C},dualnodedims::Tuple{Int,Int};dtype=Float64) where {C <: CellType}= NodePair(R,R,dualnodedims,dtype=dtype)

NodePair(R::Type{C}, nodes::Nodes{S,NX,NY,T}) where {C <: CellType, S <: CellType, NX,NY,T} = NodePair(R, (NX, NY),dtype=T)

(::Type{NodePair{R,S,NX,NY,T}})() where {R<:CellType,S<:CellType,NX,NY,T} = NodePair(R, S, (NX, NY),dtype=T)

NodePair(R, S, ::GridData{NX,NY,T};dtype=T) where {NX, NY,T} = NodePair(R, S, (NX,NY),dtype=dtype)

Base.similar(::NodePair{R,S,NX,NY,T};element_type=T) where {R<:CellType,S<:CellType,NX,NY,T} = NodePair(R, S, (NX, NY),dtype=element_type)
=#





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
