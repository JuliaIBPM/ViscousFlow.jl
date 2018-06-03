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
