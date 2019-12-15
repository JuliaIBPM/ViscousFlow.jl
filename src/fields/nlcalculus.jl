"""
    directional_derivative!(out,f,q)

Compute the directional derivative of `f` in the direction of `q`, ``q\\cdot\\nabla f``,
and put the result into `out`. Note that the result is not scaled by any grid spacing.
"""
function directional_derivative!(out::Edges{C},f::Edges{C},q::Edges{C}) where {C <: CellType}

    # some temp storage
    q_dual = Nodes(othertype(C),f)
    q_primal = Nodes(C,f)
    outx = XEdges(C,f)
    outy = YEdges(C,f)

    # compute the gradient of f
    df = EdgeGradient(C,f)
    grad!(df,f)

    # compute qy*(dfx/dy) after first transforming q.v to dual nodes and then transforming product to primal x-edges
    grid_interpolate!(outx,grid_interpolate!(q_dual,q.v) ∘ df.dudy)
    out.u .= outx

    # compute qx*(dfx/dx) after first transforming q.u to primal nodes
    # and then interpolating product to primal x-edges
    grid_interpolate!(outx,grid_interpolate!(q_primal,q.u) ∘ df.dudx)
    out.u .+= outx

    # compute qx*(dfy/dx) after first transforming qx to dual nodes and then transforming product to primal y-edges
    grid_interpolate!(outy,grid_interpolate!(q_dual,q.u) ∘ df.dvdx)
    out.v .= outy

    # compute qy*(dfy/dy) after first transforming qy to primal nodes
    # and then interpolating product to primal y-edges
    grid_interpolate!(outy,grid_interpolate!(q_primal,q.v) ∘ df.dvdy)
    out.v .+= outy

    return out

end

"""
    directional_derivative_conserve!(out,f,q)

Compute the conservative form of the directional derivative of `f` in the direction
of `q`, ``\\nabla\\cdot(qf)``, and put the result
into `out`. Note that the result is not scaled by any grid spacing. This form is only
appropriate if `q` is divergence free.
"""
function directional_derivative_conserve!(out::Edges{C},f::Edges{C},q::Edges{C}) where {C <: CellType}

    # some temp storage
    q_dual = Nodes(othertype(C),f)
    f_dual = Nodes(othertype(C),f)

    qf = q ∘ f
    qfnode = NodePair(C,othertype(C),f)

    # interpolate both fx and qx to primal nodes, multiply element by element, assign to NodePair.u
    grid_interpolate!(qfnode.u,qf.u)

    # interpolate both f.u and q.v to dual nodes, multiply element by element, assign to NodePair.v
    qfnode.v .= grid_interpolate!(q_dual,q.v) .* grid_interpolate!(f_dual,f.u)

    # take divergence and assign to out.u
    divergence!(out.u,qfnode)

    qfnode = NodePair(othertype(C),C,f)

    # interpolate both fy and qx to dual nodes, multiply element by element, assign to NodePair.u
    qfnode.u .= grid_interpolate!(q_dual,q.u) .* grid_interpolate!(f_dual,f.v)

    # interpolate both fy and qy to primal nodes, multiply element by element, assign to NodePair.v
    grid_interpolate!(qfnode.v,qf.v)

    # take divergence and assign to out.v
    divergence!(out.v,qfnode)

    return out

end

"""
    curl_cross!(out,a,b)

Compute the curl of the cross product of `a` and `b`,
``\\nabla\\times(a\\times b)``,
and put the result into `out`. Note that the result is not scaled by any grid spacing.
"""
function curl_cross!(out::Edges{C},a::Edges{C},b::Edges{C}) where {C <: CellType}

    # some temp storage
    bx_dual = Nodes(othertype(C),a)
    by_dual = Nodes(othertype(C),a)

    ax_dual = Nodes(othertype(C),a)
    ay_dual = Nodes(othertype(C),a)

    acrossb = Nodes(othertype(C),a)

    grid_interpolate!((ax_dual,ay_dual),a)
    grid_interpolate!((bx_dual,by_dual),b)

    acrossb .= ax_dual.*by_dual - ay_dual.*bx_dual
    curl!(out,acrossb)

    return out

end

"""
    convective_derivative!(out,q)

Compute the convective derivative of `q` in the form
``u\\cdot\\nabla u``
and put the result into `out`. Note that the result is not scaled by any grid spacing.
"""
convective_derivative!(out::Edges{C},q::Edges{C}) where {C<:CellType} = directional_derivative!(out,q,q)


"""
    convective_derivative_rot!(out,q)

Compute the rotational form of the convective derivative of `q` in the form
``\\frac{1}{2}\\nabla|u|^2-u\\times(\\nabla\\times u)``
and put the result into `out`. Note that the result is not scaled by any grid spacing.
"""
function convective_derivative_rot!(out::Edges{C},q::Edges{C}) where {C<:CellType}

    # some temp storage
    qx_primal = Nodes(C,q)
    qy_primal = Nodes(C,q)
    ω = curl(q)

    grid_interpolate!((qx_primal,qy_primal),q)
    grad!(out,qx_primal ∘ qx_primal + qy_primal ∘ qy_primal)
    out .*= 0.5

    q_dual = Nodes(othertype(C),q)
    qcrossω = Edges(C,q)
    grid_interpolate!(qcrossω.u,grid_interpolate!(q_dual, q.v) ∘ ω)
    grid_interpolate!(qcrossω.v,grid_interpolate!(q_dual,-q.u) ∘ ω)

    out .-= qcrossω

    return out

end
