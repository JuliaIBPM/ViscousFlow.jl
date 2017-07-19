function curl!(edges::Edges{Primal, NX, NY},
               nodes::DualNodes{NX, NY}) where {NX, NY}

    s = nodes.data
    @inbounds for y in 1:size(edges.u,2), x in 1:size(edges.u,1)
        edges.u[x,y] = s[x,y+1] - s[x,y]
    end

    @inbounds for y in 1:size(edges.v,2), x in 1:size(edges.v,1)
        edges.v[x,y] = s[x,y] - s[x+1,y]
    end
    edges
end

curl(nodes::DualNodes) = curl!(Edges(Primal, nodes), nodes)

function curl!(nodes::DualNodes{NX, NY},
               edges::Edges{Primal, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v
    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = u[x,y-1] - u[x,y] - v[x-1,y] + v[x,y]
    end
    nodes
end

function curl(edges::Edges{Primal, NX, NY}) where {NX, NY}
    curl!(DualNodes(NX, NY), edges)
end

function divergence!(nodes::DualNodes{NX, NY},
                     edges::Edges{Dual, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = - u[x-1,y-1] + u[x,y-1] - v[x-1,y-1] + v[x-1,y]
    end
    nodes
end

function divergence(edges::Edges{Dual, NX, NY}) where {NX, NY}
    divergence!(DualNodes(NX, NY), edges)
end
