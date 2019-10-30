using RecipesBase
using ColorTypes
import PlotUtils: cgrad

#using Compat
#using Compat: range

const mygreen = RGBA{Float64}(151/255,180/255,118/255,1)
const mygreen2 = RGBA{Float64}(113/255,161/255,103/255,1)
const myblue = RGBA{Float64}(74/255,144/255,226/255,1)

@recipe function plot(w::Fields.Nodes{T,NX,NY}) where {T,NX,NY}
  grid --> :none
  ratio := 1
  linewidth --> 1
  legend --> :none
  framestyle --> :frame
  levels --> range(minimum(w.data),stop=maximum(w.data),length=16)
  @series begin
    seriestype --> :contour
    transpose(w.data)
  end
end

@recipe function plot(x::AbstractArray{S,1},y::AbstractArray{S,1},w::Fields.Nodes{T,NX,NY};trim=0) where {S,T,NX,NY}
      grid --> :none
      ratio := 1
      linewidth --> 1
      legend --> :none
      framestyle --> :frame
      levels --> range(minimum(w.data),stop=maximum(w.data),length=16)
      @series begin
        seriestype --> :contour
        x[1+trim:end-trim],y[1+trim:end-trim],transpose(w.data[1+trim:end-trim,1+trim:end-trim])
      end
end

@recipe function plot(q::T) where {T <: Union{Fields.Edges,Fields.NodePair}}
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
      #levels --> range(minimum(wx.data),stop=maximum(wx.data),length=16)
      levels --> range(minimum(q.u),stop=maximum(q.u),length=16)
      #transpose(wx.data)
      transpose(q.u)
    end

    @series begin
      subplot := 2
      #levels --> range(minimum(wy.data),stop=maximum(wy.data),length=16)
      levels --> range(minimum(q.v),stop=maximum(q.v),length=16)
      #transpose(wy.data)
      transpose(q.v)
    end
end

@recipe function plot(b::Body)
    x = [b.x; b.x[1]]
    y = [b.y; b.y[1]]
    linecolor --> mygreen
    fillrange --> 0
    fillcolor --> mygreen
    ratio := 1
    legend := :none
    grid := false
    x := x
    y := y
    ()
end

function RecipesBase.RecipesBase.apply_recipe(plotattributes::Dict{Symbol, Any}, bl::BodyList)
    series_list = RecipesBase.RecipeData[]
    for b in bl
        append!(series_list, RecipesBase.RecipesBase.apply_recipe(copy(plotattributes), b) )
    end
    series_list
end
