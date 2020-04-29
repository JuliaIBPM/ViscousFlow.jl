import Base: size, ∘

#=
To do:
* Add a parameter of all GridData types that describes the type of the data array
  This will allow the data field to be of other AbstractMatrix types, like
  reshaped array
* Create constructors for all GridData types that accept an input vector
  and reshape this into an array (or arrays, via view, for Edges)
=#


macro griddata(wrapper)

  return esc(quote

    export $wrapper

    # This allows easy construction from existing GridData on the same grid.
    $wrapper(C, ::GridData{NX,NY,T};dtype=T) where {NX, NY,T <: Number} = $wrapper(C, (NX, NY),dtype=dtype )

    $wrapper(C, nx::Int, ny::Int;dtype=Float64) = $wrapper(C,(nx,ny),dtype=dtype)
    (::Type{$wrapper{C,NX,NY,T,DT}})() where {C,NX,NY,T,DT} = $wrapper(C, (NX, NY),dtype=T)

    # This constructor might be problematic? Introduced only because we have not
    # yet updated the regularization routines for the new GridData parameterization
    (::Type{$wrapper{C,NX,NY,T}})() where {C,NX,NY,T} = $wrapper(C, (NX, NY),dtype=T)

    Base.similar(::$wrapper{C,NX,NY,T,DT};element_type=T) where {C,NX,NY,T,DT} = $wrapper(C, (NX, NY),dtype=element_type)

  end)

end


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
    - `$($wrapper)`
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

    @griddata($wrapper)

    function Base.show(io::IO, u::$wrapper{C, NX, NY,T,DT}) where {C, NX, NY, T, DT}
        nodedims = "(nx = $NX, ny = $NY)"
        dims = "(nx = $(size(u,1)), ny = $(size(u,2)))"
        println(io, "$C $wrapname in a $nodedims cell grid of type $T data")
        print(io, "  Number of $C $wrapname: $dims")
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
