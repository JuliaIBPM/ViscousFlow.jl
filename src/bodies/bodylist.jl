import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect,view

export BodyList, getrange, numpts

struct BodyList
    list :: Vector{Body}
end

BodyList() = BodyList(Body[])

@propagate_inbounds getindex(A::BodyList, i::Int) = A.list[i]
@propagate_inbounds setindex!(A::BodyList, v::Body, i::Int) = A.list[i] = v
@propagate_inbounds getindex(A::BodyList, I...) = A.list[I...]
@propagate_inbounds setindex!(A::BodyList, v, I...) = A.list[I...] = v

iterate(A::BodyList) = iterate(A.list)
iterate(A::BodyList,I) = iterate(A.list,I)
size(A::BodyList) = size(A.list)
length(A::BodyList) = length(A.list)
numpts(A::BodyList) = mapreduce(b -> length(b.x),+,A)

numpts(A::Body{N}) where {N} = N

push!(bl::BodyList,b::Body) = push!(bl.list,b)

"""
    collect(bl::bodylist) -> Vector{Float64}, Vector{Float64}

Collect the inertial-space coordinates of all of the Lagrange points comprising
the bodies in body list `bl` and return each assembled set of coordinates as a vector.
"""
function collect(bl::BodyList)
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.x)
        append!(ytmp,b.y)
    end
    return xtmp,ytmp
end

collect(body::Body) = collect(BodyList([body]))

"""
    getrange(bl::BodyList,i::Int) -> Range

Return the range of indices in the global set of Lagrange point data corresponding to body `i` in body list `bl`.
"""
function getrange(bl::BodyList,i::Int)
    i <= length(bl) || error("Unavailable body")
    first = 1
    j = 1
    while j < i
        first += length(bl[j])
        j += 1
    end
    last = first+length(bl[i])-1
    return first:last
end

"""
    view(f::AbstractVector,bl::BodyList,i::Int) -> SubArray

Provide a view of the range of values in vector `f` corresponding to the Lagrange
points of the body with index `i` in body list `bl`.
"""
function Base.view(f::AbstractVector,bl::BodyList,i::Int)
    length(f) == numpts(bl) || error("Inconsistent size of data for viewing")
    return view(f,getrange(bl,i))
end

"""
    sum(f::AbstractVector,bl::BodyList,i::Int) -> Real

Compute a sum of the elements of vector `f` corresponding to body `i` in body
list `bl`.
"""
Base.sum(f::AbstractVector,bl::BodyList,i::Int) = sum(view(f,bl,i))
