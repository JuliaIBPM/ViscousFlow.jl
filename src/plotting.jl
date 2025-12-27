using Plots
pythonplot()

#=
This is for purposes of rotating the plots, from the body-fixed coordinate system in which they are solved,
to the inertial coordinate system.
=#
import PythonPlot.matplotlib as mpl
using PythonPlot.PythonCall

function get_trans_data(sol::ODESolution,sys::ILMSystem,t::Real)
    #access aux state angle, convert to degrees
    α = aux_state(sol(t))[1]*180/pi
    #access aux state x val
    x = aux_state(sol(t))[2]
    #access aux state y val
    y = aux_state(sol(t))[3]
    #perform rotation and translation
    rotationAndTranslation = mpl.transforms.Affine2D().rotate_deg(α).translate(x,y)
    
    
    trans_data = rotationAndTranslation
    
end


function pyplot_field!(ax::Py,field::Function,sol::ODESolution,sys::ILMSystem,t::Real;cmap=mpl.colormaps["RdBu"].reversed(),levels=range(-2.5,2.5,length=30))
    
    f = field(sol,sys,t)
    
    g = sys.base_cache.g
    α = aux_state(sol(t))[1]*180/pi
    x, y = coordinates(f,g)
    
    #search up matplotlib translation function
    #do rotation first and then translation
    #look at some rigid body tools (parallel axis thrm)
    #this body only has rotation
    #use the other example
    
    # we already have the plot, we want to rotate and transform the plot back to original state
    
    pts = points(sys)
    xb = [pts.u;]
    yb = [pts.v;]
    

    ax.set_aspect("equal")
    
    trans_data = get_trans_data(sol,sys,t) + ax.transData
    co = ax.contour(x,y,f',levels=levels,transform=trans_data,zorder=1,cmap=cmap)
    bo = ax.fill(xb,yb,color="black",transform=trans_data,zorder=2)
    
    co, bo
end