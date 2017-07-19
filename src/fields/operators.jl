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

function curl!(nodes::Nodes{Dual}, edges::Edges{Primal})
    u, v = edges.u, edges.v
    for y in 2:size(nodes,2)-1, x in 2:size(nodes,1)-1
        nodes.data[x,y] = u[x,y-1] - u[x,y]  + v[x,y] - v[x-1,y]
    end
    nodes
end

#function curl!(cell,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
#    @. cell[ir,jr] = -facex[ir,jr]+facex[ir,jr-1]+facey[ir,jr]-facey[ir-1,jr]
#    nothing
#end
