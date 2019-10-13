# History vectors of system solutions

import Statistics: mean
import Base: +, -, *, /

export History, HistoryType, PeriodicHistory, RegularHistory, set_first_ghost!, set_last_ghost!

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
mutable struct History{T,H <: HistoryType} <: AbstractVector{T}
  r :: UnitRange  # interior elements of vector
  vec :: Vector{T}
end

History(::Type{T}; htype::Type{H}=RegularHistory) where {T, H <: HistoryType} = History{T,htype}(1:0,T[])

History(::T;htype=RegularHistory) where {T} = History(T;htype=htype)

History(::History{T,H}) where {T, H<:HistoryType}= History(T,htype=H)

"""
    History(data::Vector[;htype=RegularHistory])

Create a history data vector and fill it with the vector `data`. The history type
can be alternatively specified as `RegularHistory` (the default) or `PeriodicHistory`
with the `htype` optional argument.
"""
History(data::Vector{T}; htype::Type{H}=RegularHistory) where {T, H <: HistoryType} = History{T,htype}(1:length(data),data)

Base.length(h::History) = length(h.vec)
Base.size(h::History) = size(h.vec)
Base.@propagate_inbounds Base.getindex(h::History, i::Int) = h.vec[i]
Base.@propagate_inbounds Base.setindex!(h::History{T}, v, i::Int) where {T} = h.vec[i] = convert(T, v)

Base.@propagate_inbounds Base.getindex(h::History{T,PeriodicHistory}, i::Int) where {T} = h.vec[mod(i-1,length(h.vec))+1]
Base.@propagate_inbounds Base.getindex(h::History{T,PeriodicHistory}, r::AbstractRange) where {T} = h.vec[mod.(r.-1,length(h.vec)).+1]

Base.push!(h::History,v...) = (h.r = h.r.start:h.r.stop+1; push!(h.vec,v...))

#=
Basic arithmetic
=#

function (-)(h_in::History)
  h = deepcopy(h_in)
  @. h.vec = -h.vec
  return h
end

# Add and subtract the same type
function (-)(h1::History{T,H},h2::History{T,H}) where {T,H}
  return History(h1.vec .- h2.vec,htype=H)
end

function (+)(h1::History{T,H},h2::History{T,H}) where {T,H}
  return History(h1.vec .+ h2.vec,htype=H)
end

# Multiply and divide by a constant
function (*)(h::History{T,H},c::Number) where {T,H}
  return History(c*h.vec,htype=H)
end


function (/)(h::History{T,H},c::Number) where {T,H}
  return History(h.vec ./ c,htype=H)
end

(*)(c::Number,h::T) where {T<:History} = *(h,c)

#=
 mean
=#
function mean!(h̄::T,h::History{T}) where {T}
    fill!(h̄,0.0)
    for hi in h
        h̄ .+= hi
    end
    h̄ /= length(h)
    return h̄
end

mean(h::History{T}) where {T} = mean!(T(),h)

#=
diff
=#
Base.diff(h::History{T,RegularHistory}) where {T} = History(diff(h.vec),htype=RegularHistory)

function Base.diff(h::History{T,PeriodicHistory}) where {T}
  r = axes(h.vec)
  r0 = ntuple(i -> i == 1 ? UnitRange(1, last(r[i])) : UnitRange(r[i]), 1)
  r1 = ntuple(i -> i == 1 ? UnitRange(2, last(r[i]) + 1) : UnitRange(r[i]), 1)

  # use unsafe_view to disable bounds checking
  return History(Base.unsafe_view(h,r1...) .- Base.unsafe_view(h,r0...), htype=PeriodicHistory)
end

Base.circshift(h::History{T,H},shift::Integer) where {T,H} = History(circshift(h.vec,shift),htype=H)

#=
set ghosts
=#
"""
    set_first_ghost!(h::History,h_pre::History)

Set the first ghost value of history `h` with the last element of history `h_pre`.
This is only valid if `h` is of type `RegularHistory`.
"""
set_first_ghost!(h::History{T,RegularHistory},h_pre::History{T}) where {T} =
    (h.r = h.r.start+1:h.r.stop+1; pushfirst!(h.vec,h_pre.vec[end]))

"""
    set_last_ghost!(h::History,h_post::History)

Set the last ghost value of history `h` with the first element of history `h_post`.
This is only valid if `h` is of type `RegularHistory`.
"""
set_last_ghost!(h::History{T,RegularHistory},h_post::History{T}) where {T} =
    push!(h.vec,h_post.vec[1])
