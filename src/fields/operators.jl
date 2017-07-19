function curl!(edges::Edges{Primal}, nodes::DualNodes)
    s = nodes.data
    for y in 1:size(edges.u,2), x in 1:size(edges.u,1)
        edges.u[x,y] = s[x,y+1] - s[x,y]
    end

    for y in 1:size(edges.v,2), x in 1:size(edges.v,1)
        edges.v[x,y] = s[x,y] - s[x+1,y]
    end
    edges
end

curl(nodes::DualNodes) = curl!(Edges(Primal, nodes), nodes)

function curl!(nodes::DualNodes, edges::Edges{Primal})
    u, v = edges.u, edges.v
    for y in 2:size(nodes,2)-1, x in 2:size(nodes,1)-1
        nodes[x,y] = u[x,y-1] - u[x,y] - v[x-1,y] + v[x,y]
    end
    nodes
end

curl(edges::Edges{Primal}) = curl!(DualNodes(edges.nodedims), edges)
