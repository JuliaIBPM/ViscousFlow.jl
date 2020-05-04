import Base: size, ∘

"""
    @scalarfield(wrapper,primaldn,dualdn)

Given a name `wrapper`, generate a scalar grid data type that wraps an
array of grid data. The macro automatically generates several constructors
and functions for the data type. The tuples `primaldn` and `dualdn` specify
the different number of indices in each direction compared to the reference
grid type (dual nodes). For example, if `primaldn` is (-1,0), then this corresponds
to one fewer point in the x direction and the same number of points in the y
direction compared to dual nodes. Each 0 entry in these tuples will be interpreted
as placing the data half a cell spacing from the boundary of the domain in
the physical grid.
"""
macro scalarfield(wrapper,primaldn,dualdn)

  wrapname = lowercase(string(wrapper))
  indexfcn = Symbol(wrapname[1:end-1],"_inds")

  #pdn = 0 .- eval(primaldn)
  #ddn = 0 .- eval(dualdn)
  #pshift = 0.5.*(1 .- pdn)
  #dshift = 0.5.*(1 .- ddn)

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

    # functions set the number of indices in each direction, based on the
    # number of dual node dimensions
    export $indexfcn

    $indexfcn(::Type{Dual},   dualnodedims) = dualnodedims .+ $dualdn
    $indexfcn(::Type{Primal}, dualnodedims) = dualnodedims .+ $primaldn

    Base.size(::Type{<:$wrapper{C,NX,NY}}) where {C,NX,NY} = $indexfcn(C,(NX,NY))

    # Provide the correction to be applied to the reference grid's number
    # of indices in each direction
    indexshift(::$wrapper{Dual}) = $dualdn
    indexshift(::$wrapper{Primal}) = $primaldn
    indexshift(::Type{T}) where {T<:$wrapper{Dual,NX,NY,T,DT} where {NX,NY,T,DT}} = $dualdn
    indexshift(::Type{T}) where {T<:$wrapper{Primal,NX,NY,T,DT} where {NX,NY,T,DT}} = $primaldn

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
