## Field data macros

export show_scalarlist, show_vectorlist, show_tensorlist

_negate(a::NTuple{M,Int64}) where {M} = map(x -> -x,a)
_cellshift(a::NTuple{M,Int64}) where {M} = map(x -> 0.5*(1-abs(x)),a)

# for Julia < 1.1
_fieldtypes(T::Type) = ntuple(i -> fieldtype(T, i), fieldcount(T))

"""
    @griddata(wrapper,nctypes)

Create a basic set of constructors and functions for grid data wrapper type
`wrapper`. The argument `nctypes` is the number of cell types used
to parameterize this grid data type.
"""
macro griddata(wrapper, nctypes)

  ctype = Symbol[]
  for i in 1:eval(nctypes)
    push!(ctype,Symbol("C",string(i)))
  end

  return esc(quote

    export $wrapper

    celltype(::$wrapper{C}) where {C<:CellType} = C
    celltype(::Type{<:$wrapper{$(ctype...),NX,NY}}) where {$(ctype...),NX,NY} = $(ctype...)

    gridsize(::Type{<:$wrapper{$(ctype...),NX,NY}}) where {$(ctype...),NX,NY} = NX, NY

    griddatatype(::$wrapper) = $wrapper
    griddatatype(::Type{T}) where {T<:$wrapper{$(ctype...),NX,NY,T,DT} where {$(ctype...),NX,NY,T,DT}} = $wrapper

    # This allows easy construction from existing GridData on the same grid.
    $wrapper(C, ::GridData{NX,NY,T};dtype=T) where {NX, NY,T <: Number} = $wrapper(C, (NX, NY),dtype=dtype )

    $wrapper(C, nx::Int, ny::Int;dtype=Float64) = $wrapper(C,(nx,ny),dtype=dtype)
    (::Type{$wrapper{$(ctype...),NX,NY,T,DT}})() where {$(ctype...),NX,NY,T,DT} =
                  $wrapper($(ctype[1]), (NX, NY),dtype=T)

    # This constructor might be problematic? Introduced only because we have not
    # yet updated the regularization routines for the new GridData parameterization
    (::Type{$wrapper{$(ctype...),NX,NY,T}})() where {$(ctype...),NX,NY,T} =
                  $wrapper($(ctype[1]), (NX, NY),dtype=T)

    Base.similar(::$wrapper{$(ctype...),NX,NY,T,DT};element_type=T) where {$(ctype...),NX,NY,T,DT} =
                  $wrapper($(ctype[1]), (NX, NY),dtype=element_type)

    Base.parent(A::$wrapper) = A.data
    Base.size(A::$wrapper) = size(A.data)
    Base.parentindices(A::$wrapper) = parentindices(A.data)

    @propagate_inbounds Base.getindex(A::$wrapper,i::Int) = getindex(A.data,i)
    @propagate_inbounds Base.setindex!(A::$wrapper, v, i::Int) = setindex!(A.data,convert(eltype(A.data),v),i)

  end)

end

"""
    @generate_scalarlist(list)

Process a list of scalar grid types and generate an expanded form of the list,
complete with cell shifts. Used for constructing the regularization functions
and coordinate functions.
"""
macro generate_scalarlist(list)
    return quote
        newlist = []
        for (i,l) in enumerate($(esc(list)))
            wrapper, primaldn, dualdn = l

            # negative index shifts become positive values in this list
            pdn = _negate(primaldn)
            ddn = _negate(dualdn)

            # grid cell shifts: dn = 0 implies shift by half-cell, dn = 1 implies no shift
            pshift = _cellshift(pdn)
            dshift = _cellshift(ddn)
            push!(newlist,(Symbol(wrapper),:Primal,pdn...,pshift...))
            push!(newlist,(Symbol(wrapper),:Dual,ddn...,dshift...))
        end
        newlist
    end
end

show_scalarlist() = @generate_scalarlist SCALARLIST

"""
    @generate_vectorlist(list)

Process a list of vector grid types and generate an expanded form of the list,
complete with cell shifts. Used for constructing the regularization functions
and coordinate functions.
"""
macro generate_collectionlist(list)
    return quote
        newlist = []
        for (i,l) in enumerate($(esc(list)))
            wrapper, celltypes = l
            gtype = eval(wrapper){eval.(celltypes)...}

            row = (l...,)
            for ft in unique(_fieldtypes(gtype))
                if ft <: GridData
                    dn = _negate(indexshift(ft))
                    row = (row...,dn...)
                end
            end
            for ft in unique(_fieldtypes(gtype))
                if ft <: GridData
                    dn = _negate(indexshift(ft))
                    cshift = _cellshift(dn)
                    row = (row...,cshift...)
                end
            end
            push!(newlist,row)
        end
        newlist
    end
end

show_vectorlist() = @generate_collectionlist VECTORLIST

show_tensorlist() = @generate_collectionlist TENSORLIST


"""
    @wrapparay(wrapper,field,N)

Basic macro to develop any AbstractArray data type into a proper wrapper, with
indexing and other array-type operations.
"""
macro wraparray(wrapper, field, N)
    S = eval(wrapper)
    @assert S <: AbstractArray "Wrapped type must be a subtype of AbstractArray"

    quote
        Base.parent(A::$wrapper) = A.$field
        Base.size(A::$wrapper) = size(A.$field)
        parentindices(A::$wrapper) = parentindices(A.$field)

        if $N > 1
          function Base.show(io::IO, m::MIME"text/plain", A::$wrapper)
            println(io, "$(typeof(A)) data")
            println(io, "Printing in grid orientation (lower left is (1,1))")
            show(io,m, reverse(transpose(A.$field),dims=1))
          end
          #function Base.summary(io::IO, A::$wrapper)
          #  println(io, "$(typeof(A)) data")
          #  print(io, "Printing in grid orientation (lower left is (1,1))")
          #end
        end

        @propagate_inbounds Base.getindex(A::$wrapper, i::Int) = A.$field[i]
        @propagate_inbounds Base.setindex!(A::$wrapper, v, i::Int) = A.$field[i] = convert(eltype(A.$field), v)
        if $N > 1
          @propagate_inbounds Base.getindex(A::$wrapper, I::Vararg{Int, $N}) = A.$field[I...]
          @propagate_inbounds Base.setindex!(A::$wrapper, v, I::Vararg{Int, $N}) = A.$field[I...] = convert(eltype(A.$field), v)
        end
    end
end
