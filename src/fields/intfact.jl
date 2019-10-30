# Integrating factor

intfact(x, a) = exp(-2a)besseli(x,2a)

intfact(x, y,a) = exp(-4a)besseli(x,2a)besseli(y,2a)

intfact(x, y, z, a) = exp(-6a)besseli(x,2a)besseli(y,2a)besseli(z,2a)


"""
    plan_intfact(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])

Constructor to set up an operator for evaluating the integrating factor with
real-valued parameter `a`. This can then be applied with the `*` operation on
data of the appropriate size.

The `dims` argument can be replaced with `w::Nodes` to specify the size of the
domain.

# Example

```jldoctest
julia> w = Nodes(Dual,(6,6));

julia> w[4,4] = 1.0;

julia> E = plan_intfact(1.0,(6,6))
Integrating factor with parameter 1.0 on a (nx = 6, ny = 6) grid

julia> E*w
Nodes{Dual,6,6,Float64} data
Printing in grid orientation (lower left is (1,1))
6×6 Array{Float64,2}:
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00888233   0.028765    0.0664124   0.0951774   0.0664124   0.028765
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.000828935  0.00268447  0.00619787  0.00888233  0.00619787  0.00268447
```
"""
function plan_intfact end

"""
    plan_intfact!(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])

Same as [`plan_intfact`](@ref), but the resulting operator performs an in-place
operation on data.
"""
function plan_intfact! end


struct IntFact{NX, NY, a, inplace}
    conv::Union{CircularConvolution{NX, NY},Nothing}
end

for (lf,inplace) in ((:plan_intfact,false),
                     (:plan_intfact!,true))

    @eval function $lf(a::Real,dims::Tuple{Int,Int};fftw_flags = FFTW.ESTIMATE)
        NX, NY = dims

        if a == 0
          return IntFact{NX, NY, 0.0, $inplace}(nothing)
        end

        #qtab = [intfact(x, y, a) for x in 0:NX-1, y in 0:NY-1]
        Nmax = 0
        while abs(intfact(Nmax,0,a)) > eps(Float64)
          Nmax += 1
        end
        qtab = [max(x,y) <= Nmax ? intfact(x, y, a) : 0.0 for x in 0:NX-1, y in 0:NY-1]
        #IntFact{NX, NY, a, $inplace}(Nullable(CircularConvolution(qtab, fftw_flags)))
        IntFact{NX, NY, a, $inplace}(CircularConvolution(qtab, fftw_flags))
      end

      @eval $lf(a::Real,nodes::Nodes{T,NX,NY}; fftw_flags = FFTW.ESTIMATE) where {T<:CellType,NX,NY} =
          $lf(a,node_inds(T,(NX,NY)), fftw_flags = fftw_flags)


end


function Base.show(io::IO, E::IntFact{NX, NY, a, inplace}) where {NX, NY, a, inplace}
    nodedims = "(nx = $NX, ny = $NY)"
    isinplace = inplace ? "In-place integrating factor" : "Integrating factor"
    print(io, "$isinplace with parameter $a on a $nodedims grid")
end

function mul!(out::Nodes{T,NX, NY},
                   E::IntFact{MX, MY, a, inplace},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, a, inplace}

    mul!(out.data, E.conv, s.data)
    out
end

function mul!(out::Nodes{T,NX, NY},
                   E::IntFact{MX, MY, 0.0, inplace},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, inplace}
    out .= deepcopy(s)
end

*(E::IntFact{MX,MY,a,false},s::Nodes{T,NX,NY}) where {MX,MY,a,T <: CellType, NX,NY} =
  mul!(Nodes(T,s), E, s)

*(E::IntFact{MX,MY,a,true},s::Nodes{T,NX,NY}) where {MX,MY,a,T <: CellType, NX,NY} =
    mul!(s, E, deepcopy(s))
