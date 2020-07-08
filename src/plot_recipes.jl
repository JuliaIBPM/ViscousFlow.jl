using RecipesBase
using ColorTypes
using LaTeXStrings
import PlotUtils: cgrad


@recipe f(w::T,sys::NavierStokes) where {T<:GridData} = w, sys.grid
