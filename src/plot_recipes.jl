using RecipesBase
using ColorTypes
import PlotUtils: cgrad

@recipe function plot(w::Fields.Nodes{T,NX,NY}) where {T,NX,NY}
    seriestype --> :contour
      grid --> :none
      ratio := 1
      linewidth --> 1
      legend --> :none
      framestyle --> :frame
      levels --> linspace(minimum(w.data),maximum(w.data),16)
      transpose(w.data)
end

@recipe function plot(q::Fields.Edges{T,NX,NY}) where {T,NX,NY}
    #wx = Fields.Nodes(Dual,(NX,NY))
    #wy = Fields.Nodes(Dual,(NX,NY))
    #shift!((wx,wy),q)
    layout := (1,2)
    seriestype --> :contour
    grid --> :none
    ratio := 1
    linewidth --> 1
    legend --> :none
    framestyle --> :frame
    @series begin
      subplot := 1
      #levels --> linspace(minimum(wx.data),maximum(wx.data),16)
      levels --> linspace(minimum(q.u),maximum(q.u),16)
      #transpose(wx.data)
      transpose(q.u)
    end

    @series begin
      subplot := 2
      #levels --> linspace(minimum(wy.data),maximum(wy.data),16)
      levels --> linspace(minimum(q.v),maximum(q.v),16)
      #transpose(wy.data)
      transpose(q.v)
    end
end
