import Base: size, ∘


macro scalarfield(wrapper)

  wrapname = lowercase(string(wrapper))
  indexfcn = Symbol(wrapname[1:end-1],"_inds")

  return esc(quote

    # The data type

    @doc """
        $($wrapper)

    `$($wrapper)` is a wrapper for scalar-valued data that lie at the centers of either dual cells or
    primary cells. A `$($wrapper)` type can be accessed by indexing like any other array,
    and allows the use of [`size`](@ref), [`similar`](@ref), [`zero`](@ref).

    # Constructors
    - `$($wrapper)(C,dims)` creates a field of zeros in cells of type `C` (where `C` is
      either `Dual` or `Primal`), on a grid of dimensions `dims` (a tuple). Note that `dims`
      represent the number of dual cells on the grid, even if `C` is `Primal`.
    - `$($wrapper)(C,w)` performs the same construction, but uses existing field data `w`
      of `GridData` type to determine the size of the grid.
    -  Adding the `dtype=` keyword allows the data type of the field data to be
      changed. The default is `Float64`, but can be changed to, e.g., `ComplexF64`
    """
    struct $wrapper{C <: CellType, NX, NY, T <: Number, DT <: AbstractMatrix} <: ScalarGridData{NX,NY,T}
      data::DT
      $wrapper{C,NX,NY,T,DT}(data::AbstractMatrix) where {C<: CellType,NX,NY,T<:Number,DT} =
            new{C,NX,NY,T,typeof(data)}(data)
      $wrapper{C,NX,NY,T,DT}(data::AbstractVector) where {C<: CellType,NX,NY,T<:Number,DT} =
            (u = reshape(data,$indexfcn(C,(NX,NY))); new{C,NX,NY,T,typeof(u)}(u))
    end

    # Constructors

    function $wrapper(::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
        dims = $indexfcn(C, dualnodedims)
        $wrapper{C, dualnodedims...,dtype,typeof(zeros(dtype,dims))}(zeros(dtype,dims))
    end

    # Endow this scalar grid data type with some basic constructors and functions
    @griddata($wrapper,1)

    # Multidimensional indexing for scalar grid data
    @propagate_inbounds Base.getindex(A::$wrapper, I::Vararg{Int,2}) = A.data[I...]
    @propagate_inbounds Base.setindex!(A::$wrapper, v, I::Vararg{Int,2}) = A.data[I...] = convert(eltype(A.data), v)

    function Base.show(io::IO, u::$wrapper{C, NX, NY,T,DT}) where {C, NX, NY, T, DT}
        nodedims = "(nx = $NX, ny = $NY)"
        dims = "(nx = $(size(u,1)), ny = $(size(u,2)))"
        println(io, "$C $($wrapname) in a $nodedims cell grid of type $T data")
        print(io, "  Number of $C $($wrapname): $dims")
    end

    function Base.show(io::IO, m::MIME"text/plain", A::$wrapper)
      println(io, "$(typeof(A)) data")
      println(io, "Printing in grid orientation (lower left is (1,1))")
      show(io,m, reverse(transpose(A.data),dims=1))
    end


  end)

end



# [X]_inds functions set the number of indices in each direction
# Based on number of dual nodes, return the number of elements of this type

@scalarfield Nodes
node_inds(::Type{Dual},   dualnodedims) = dualnodedims[1], dualnodedims[2]
node_inds(::Type{Primal}, dualnodedims) = dualnodedims[1]-1, dualnodedims[2]-1

@scalarfield XEdges
xedge_inds(::Type{Dual}, dualnodedims) = dualnodedims[1]-1, dualnodedims[2]
xedge_inds(::Type{Primal}, dualnodedims) = dualnodedims[1], dualnodedims[2]-1

@scalarfield YEdges
yedge_inds(::Type{Dual}, dualnodedims) = dualnodedims[1], dualnodedims[2]-1
yedge_inds(::Type{Primal}, dualnodedims) = dualnodedims[1]-1, dualnodedims[2]
