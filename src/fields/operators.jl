function curl!(edges::Edges{Primal}, nodes::Nodes{Dual})
    s = nodes.data
    for y in 2:size(edges.u,2), x in 2:size(edges.u,1)
        edges.u[x,y] = s[x,y+1] - s[x,y]
    end

    for y in 2:size(edges.v,2), x in 2:size(edges.v,1)
        edges.v[x,y] = s[x,y] - s[x+1,y]
    end
    edges
end

function curl(nodes::Nodes{Dual})
    edges = Edges(Primal, nodes.celldims, nodes.ghostlayers)
    curl!(edges, nodes)
end

function curl!(nodes::Nodes{Dual}, edges::Edges{Primal})
    u, v = edges.u, edges.v
    for y in 2:size(nodes,2)-1, x in 2:size(nodes,1)-1
        nodes.data[x,y] = u[x,y-1] - u[x,y]  + v[x,y] - v[x-1,y]
    end
    nodes
end

function curl(edges::Edges{Primal})
    nodes = Nodes(Dual, edges.celldims, edges.ghostlayers)
    curl!(nodes, edges)
end
