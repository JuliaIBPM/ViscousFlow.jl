# History vectors of system solutions

import Base: length

abstract type HistoryType end

abstract type PeriodicHistory <: HistoryType end

abstract type RegularHistory <: HistoryType end

"""
    History(datatype,[;htype::HistoryType = RegularHistory])

Create a blank history data vector with entries of type `datatype`. Alternatively,
one can pass an example instance of the type of entry. An optional argument `htype`
specifies the type of history vector. By default, this is `RegularHistory`, but
can be alternatively set to `PeriodicHistory`. In the latter case, if the
history vector has length `n`, then it will be assumed that the `n+1` entry
is identical to the `1` entry.

It is important to note that `datatype` must be outfitted with basic operations:
`+`, `-`, scalar multiplication, and `fill!`.
"""
struct History{T,H <: HistoryType}
  vec :: Vector{T}
end

function History(::Type{T}; htype::Type{H}=RegularHistory) where {T, H <: HistoryType}
    println()
    return History{T,htype}(T[])
end

History(::T;htype=RegularHistory) where {T} = History(T;htype=htype)

length(h::History) = length(h.vec)

import Statistics: mean

function mean!(q̄::T,q::Vector{T}) where {T}
    fill!(q̄,0.0)
    for qi in q
        q̄ .+= qi
    end
    return q̄/length(q)
end

mean(q::Vector{T}) where {T <: GridData} = mean!(T(),q)
