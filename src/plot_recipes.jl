using RecipesBase
using ColorTypes
using LaTeXStrings
import PlotUtils: cgrad

#using Compat
#using Compat: range

const mygreen = RGBA{Float64}(151/255,180/255,118/255,1)
const mygreen2 = RGBA{Float64}(113/255,161/255,103/255,1)
const myblue = RGBA{Float64}(74/255,144/255,226/255,1)


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

@recipe function f(m::RigidBodyMotion;tmax=10)

    t = 0.0:0.01:tmax
    ux = map(ti -> real(m(ti)[2]),t)
    uy = map(ti -> imag(m(ti)[2]),t)
    adot = map(ti -> m(ti)[5],t)
    xlim --> (0,tmax)
  layout := 3
  grid --> :none
  linewidth --> 1
  legend --> :none
  framestyle --> :frame
  xlabel --> L"t"
  ulim = min(minimum(ux),minimum(uy)),max(maximum(ux),maximum(uy))
  alim = extrema(adot)
  @series begin
      subplot := 1
      ylim --> ulim
      ylabel --> L"u"
      t, ux
    end
  @series begin
      subplot := 2
      ylim --> ulim
      ylabel --> L"v"
      t, uy
    end
   @series begin
      subplot := 3
      ylim --> alim
      ylabel --> L"\dot{\alpha}"
      t, adot
    end
end

function RecipesBase.RecipesBase.apply_recipe(plotattributes::Dict{Symbol, Any}, bl::BodyList)
    series_list = RecipesBase.RecipeData[]
    for b in bl
        append!(series_list, RecipesBase.RecipesBase.apply_recipe(copy(plotattributes), b) )
    end
    series_list
end
