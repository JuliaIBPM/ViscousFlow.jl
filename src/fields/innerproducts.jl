using LinearAlgebra
import LinearAlgebra:dot, norm

"""
    dot(p1::Nodes{Dual},p2::Nodes{Dual}) -> Real

Computes the inner product between two sets of dual node data on the same grid.
"""
function dot(p1::Nodes{Dual,NX,NY},p2::Nodes{Dual,NX,NY}) where {NX,NY}
  # remember that sizes NX and NY include the ghost cells
  dims = node_inds(Dual,(NX,NY))
  return dot(p1[2:dims[1]-1,2:dims[2]-1],p2[2:dims[1]-1,2:dims[2]-1])/((NX-2)*(NY-2))
end

"""
    dot(p1::Nodes{Primal},p2::Nodes{Primal}) -> Real

Computes the inner product between two sets of primal node data on the same grid.
"""
function dot(p1::Nodes{Primal,NX,NY},p2::Nodes{Primal,NX,NY}) where {NX,NY}

  dims = node_inds(Primal,(NX,NY))

  # interior
  tmp = dot(p1[2:dims[1]-1,2:dims[2]-1],p2[2:dims[1]-1,2:dims[2]-1])

  # boundaries
  tmp += 0.5*dot(p1[1,2:dims[2]-1],      p2[1,2:dims[2]-1])
  tmp += 0.5*dot(p1[dims[1],2:dims[2]-1],p2[dims[1],2:dims[2]-1])
  tmp += 0.5*dot(p1[2:dims[1]-1,1],      p2[2:dims[1]-1,1])
  tmp += 0.5*dot(p1[2:dims[1]-1,dims[2]],p2[2:dims[1]-1,dims[2]])

  # corners
  tmp += 0.25*(p1[1,1]*p2[1,1]             + p1[dims[1],1]*p2[dims[1],1] +
               p1[1,dims[2]]*p2[1,dims[2]] + p1[dims[1],dims[2]]*p2[dims[1],dims[2]])

  return tmp/((NX-2)*(NY-2))
end

"""
    dot(p1::Edges{Dual},p2::Edges{Dual}) -> Real

Computes the inner product between two sets of dual edge data on the same grid.
"""
function dot(p1::Edges{Dual,NX,NY},p2::Edges{Dual,NX,NY}) where {NX,NY}

  udims, vdims = edge_inds(Dual,(NX,NY))

  # interior
  tmp = dot(p1.u[2:udims[1]-1,2:udims[2]-1],p2.u[2:udims[1]-1,2:udims[2]-1]) +
        dot(p1.v[2:vdims[1]-1,2:vdims[2]-1],p2.v[2:vdims[1]-1,2:vdims[2]-1])

  # boundaries
  tmp += 0.5*dot(p1.u[1,       2:udims[2]-1],p2.u[1,       2:udims[2]-1])
  tmp += 0.5*dot(p1.u[udims[1],2:udims[2]-1],p2.u[udims[1],2:udims[2]-1])
  tmp += 0.5*dot(p1.v[2:vdims[1]-1,1],       p2.v[2:vdims[1]-1,1])
  tmp += 0.5*dot(p1.v[2:vdims[1]-1,vdims[2]],p2.v[2:vdims[1]-1,vdims[2]])

  return tmp/((NX-2)*(NY-2))
end

"""
    dot(p1::Edges{Primal},p2::Edges{Primal}) -> Real

Computes the inner product between two sets of primal edge data on the same grid.
"""
function dot(p1::Edges{Primal,NX,NY},p2::Edges{Primal,NX,NY}) where {NX,NY}

  udims, vdims = edge_inds(Primal,(NX,NY))

  # interior
  tmp = dot(p1.u[2:udims[1]-1,2:udims[2]-1],p2.u[2:udims[1]-1,2:udims[2]-1]) +
        dot(p1.v[2:vdims[1]-1,2:vdims[2]-1],p2.v[2:vdims[1]-1,2:vdims[2]-1])

  # boundaries
  tmp += 0.5*dot(p1.u[2:udims[1]-1,1],       p2.u[2:udims[1]-1,1])
  tmp += 0.5*dot(p1.u[2:udims[1]-1,udims[2]],p2.u[2:udims[1]-1,udims[2]])
  tmp += 0.5*dot(p1.v[1,       2:vdims[2]-1],p2.v[1,       2:vdims[2]-1])
  tmp += 0.5*dot(p1.v[vdims[1],2:vdims[2]-1],p2.v[vdims[1],2:vdims[2]-1])

  return tmp/((NX-2)*(NY-2))
end

"""
    norm(p::GridData) -> Real

Computes the L2 norm of data on a grid.
"""
norm(p::GridData) = sqrt(dot(p,p))

# This function computes an integral by just taking the inner product with
# another set of cell data uniformly equal to 1
"""
    integrate(p::Nodes) -> Real

Computes a numerical quadrature of node data.
"""
function integrate(p::Nodes{C,NX,NY}) where {C,NX,NY}
  p2 = zero(p)
  fill!(p2.data,1) # fill it with ones
  return dot(p,p2)
end
