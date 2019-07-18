import Base: -, +, *, /, ∘

### On scalar grid data ####

# Set it to negative of itself
function (-)(p_in::ScalarGridData)
  p = deepcopy(p_in)
  p.data .= -p.data
  return p
end

function (-)(p1::T,p2::T) where {T <: ScalarGridData}
   return T(p1.data .- p2.data)
 end

function (+)(p1::T,p2::T) where {T <: ScalarGridData}
  return T(p1.data .+ p2.data)
end

# Multiply and divide by a constant
function (*)(p::T,c::Number) where {T<:ScalarGridData}
  return T(c*p.data)
end

function (/)(p::T,c::Number) where {T<:ScalarGridData}
  return T(p.data ./ c)
end

function product!(out::Nodes{T, NX, NY},
                  p::Nodes{T, NX, NY},
                  q::Nodes{T, NX, NY}) where {T, NX, NY}

    inds = node_inds(T, (NX, NY))
    @inbounds for y in 1:inds[2], x in 1:inds[1]
        out[x,y] = p[x,y] * q[x,y]
    end
    out
end

function product(p::Nodes{T, NX, NY}, q::Nodes{T, NX, NY}) where {T, NX, NY}
    product!(Nodes(T, (NX, NY)), p, q)
end

function (∘)(p::Nodes{T, NX, NY}, q::Nodes) where {T, NX, NY}
    product!(Nodes(T, (NX, NY)), p, q)
end

(*)(c::Number,p::T) where {T<:GridData} = *(p,c)


### On vector grid data ####

function (-)(p_in::VectorGridData)
  p = deepcopy(p_in)
  p.u .= -p.u
  p.v .= -p.v
  return p
end

function (-)(p1::T,p2::T) where {T<:VectorGridData}
  return T(p1.u - p2.u, p1.v - p2.v)
end

function (+)(p1::T,p2::T) where {T<:VectorGridData}
  return T(p1.u + p2.u, p1.v + p2.v)
end

function (*)(p::T,c::Number) where {T<:VectorGridData}
  return T(c*p.u,c*p.v)
end

function (/)(p::T,c::Number) where {T<:VectorGridData}
  return T(p.u / c, p.v / c)
end


"""
    product!(out::Edges/Nodes,p::Edges/Nodes,q::Edges/Nodes)

Compute the Hadamard (i.e. element by element) product of edge or nodal
(primal or dual) data `p` and `q` and return the result in `out`.

# Example

```jldoctest
julia> q = Edges(Dual,(8,6));

julia> out = p = deepcopy(q);

julia> q.u[3,2] = 0.3;

julia> p.u[3,2] = 0.2;

julia> product!(out,p,q)
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.06  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function product!(out::Edges{T, NX, NY},
                  p::Edges{T, NX, NY},
                  q::Edges{T, NX, NY}) where {T, NX, NY}

    uinds, vinds = edge_inds(T, (NX, NY))
    @inbounds for y in 1:uinds[2], x in 1:uinds[1]
        out.u[x,y] = p.u[x,y] * q.u[x,y]
    end

    @inbounds for y in 1:vinds[2], x in 1:vinds[1]
        out.v[x,y] = p.v[x,y] * q.v[x,y]
    end
    out
end

"""
    product(p::Edges/Nodes,q::Edges/Nodes) --> Edges/Nodes

Compute the Hadamard product of edge or nodal (primal or dual) data `p` and `q` and return
the result. This operation can also be carried out with the `∘` operator:

# Example

```jldoctest
julia> q = Edges(Dual,(8,6));

julia> p = deepcopy(q);

julia> q.u[3,2] = 0.3;

julia> p.u[3,2] = 0.2;

julia> p∘q
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.06  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function product(p::Edges{T, NX, NY}, q::Edges{T, NX, NY}) where {T, NX, NY}
    product!(Edges(T, (NX, NY)), p, q)
end

function (∘)(p::Edges{T, NX, NY}, q::Edges) where {T, NX, NY}
    product!(Edges(T, (NX, NY)), p, q)
end

#### ON ALL TYPES ####
