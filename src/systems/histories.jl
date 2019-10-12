# History vectors of system solutions

import Statistics: mean

export History, PeriodicHistory, RegularHistory

abstract type HistoryType end

abstract type PeriodicHistory <: HistoryType end

abstract type RegularHistory <: HistoryType end

"""
    History(datatype,[;htype::HistoryType = RegularHistory])

Create an empty history data vector with entries of type `datatype`. Alternatively,
one can pass an example instance of the type of entry. An optional argument `htype`
specifies the type of history vector. By default, this is `RegularHistory`, but
can be alternatively set to `PeriodicHistory`. In the latter case, if the
history vector has length `n`, then it will be assumed that the `n+1` entry
is identical to the `1` entry.

It is important to note that, in order to use the routines for `History` types,
then the element type `datatype` must be outfitted with basic operations:
`+`, `-`, scalar multiplication, and `fill!`.

Another constructor is `History(h::History)`, which creates an empty instance of
a history of the same type as `h`.
"""
struct History{T,H <: HistoryType} <: AbstractVector{T}
  vec :: Vector{T}
end

History(::Type{T}; htype::Type{H}=RegularHistory) where {T, H <: HistoryType} = History{T,htype}(T[])

History(::T;htype=RegularHistory) where {T} = History(T;htype=htype)

History(::History{T,H}) where {T, H<:HistoryType}= History(T,htype=H)

"""
    History(data::Vector[;htype=RegularHistory])

Create a history data vector and fill it with the vector `data`. The history type
can be alternatively specified as `RegularHistory` (the default) or `PeriodicHistory`
with the `htype` optional argument.
"""
History(data::Vector{T}; htype::Type{H}=RegularHistory) where {T, H <: HistoryType} = History{T,htype}(data)

Base.length(h::History) = length(h.vec)
Base.size(h::History) = size(h.vec)
@propagate_inbounds Base.getindex(h::History, i::Int) = h.vec[i]
@propagate_inbounds Base.setindex!(h::History{T}, v, i::Int) = h.vec[i] = convert(T, v)


function mean!(h̄::T,h::History{T}) where {T}
    fill!(h̄,0.0)
    for hi in h
        h̄ .+= hi
    end
    return h̄/length(h)
end

mean(h::History{T}) where {T} = mean!(T(),h)

Base.diff(h::History{T,RegularHistory}) where {T} = History(diff(h.vec),htype=RegularHistory)
