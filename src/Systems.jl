module Systems

import Whirl2d
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.ddf

mutable struct Domain
    xmin::Vector{Float64}	
    xmax::Vector{Float64}	

    grid::AbstractArray{Grids.DualPatch}
    ngrid::Int

    body::AbstractArray{Bodies.Body}
    nbody::Int
    nbodypts::Int
    firstbpt::Vector{Int}

    ddf_fcn

    """
    E^T operator, acting on body Lagrange points and returning data at
    cell edges
    """
    Extrans::SparseMatrixCSC{Float64,Int}
    Eytrans::SparseMatrixCSC{Float64,Int}

    """
    (EC)^T operator, acting on body Lagrange points and returning dual grid cell data
    """
    ExCtrans::SparseMatrixCSC{Float64,Int}
    EyCtrans::SparseMatrixCSC{Float64,Int}

end

function Domain()
    xmin = zeros(Whirl2d.ndim)
    xmax = ones(Whirl2d.ndim)	 
    ngrid = nbody = 0
    nbodypts = 0
    firstbpt = []
    Extrans = spzeros(0)
    Eytrans = spzeros(0)
    ExCtrans = spzeros(0)
    EyCtrans = spzeros(0)
    
    ddf_fcn = Whirl2d.ddf.ddf_roma

    Domain(xmin,xmax,[],ngrid,[],nbody,nbodypts,firstbpt,
           ddf_fcn,Extrans,Eytrans,ExCtrans,EyCtrans)

end

Domain(g::Grids.DualPatch) = add_grid(Domain(),g)
Domain(b::Bodies.Body) = add_body(Domain(),b)
Domain(xmin,xmax) = set_dims(Domain(),xmin,xmax)


function set_dims!(dom::Domain,
		  xmin::Vector{Float64},xmax::Vector{Float64})
    dom.xmin = xmin
    dom.xmax = xmax
end

function set_dims(dom::Domain,
		  xmin::Vector{Float64},xmax::Vector{Float64})
    set_dims!(dom,xmin,xmax)
    dom
end

function add_grid!(dom::Domain,g::Grids.DualPatch)
    dom.ngrid += 1
    resize!(dom.grid,dom.ngrid)
    dom.grid[dom.ngrid] = g

    # resize the domain dimensions to allow the grid to fit
    xmin = [min(dom.xmin[j],g.xmin[j]) for j = 1:Whirl2d.ndim]
    xmax = [max(dom.xmax[j],g.xmax[j]) for j = 1:Whirl2d.ndim]
    set_dims!(dom,xmin,xmax)

end

function add_grid(dom::Domain,g::Grids.DualPatch)
    add_grid!(dom,g)
    dom
end

function add_grid(dom::Domain,Δx::Float64)
    # create a domain-filling patch
    N = ceil.(Int,(dom.xmax-dom.xmin)/Δx)
    g = Grids.DualPatch(N,Δx,dom.xmin)
    add_grid!(dom,g)
    dom
end

function add_body!(dom::Domain,b::Bodies.Body)
    dom.nbody += 1
    resize!(dom.body,dom.nbody)
    dom.body[dom.nbody] = b

    resize!(dom.firstbpt,dom.nbody)
    dom.firstbpt[dom.nbody] = dom.nbodypts+1
    dom.nbodypts += b.N
    
    buff = 0.0

    # resize the domain dimensions to allow the body to fit, with a buffer
    bmin, bmax = Bodies.dims(b)
    xmin = [min(dom.xmin[j],bmin[j]-buff) for j = 1:Whirl2d.ndim]
    xmax = [max(dom.xmax[j],bmax[j]+buff) for j = 1:Whirl2d.ndim]
    set_dims!(dom,xmin,xmax)
end	

function add_body(dom::Domain,b::Bodies.Body)
    add_body!(dom,b)
    dom
end

# Body-to-grid operations
function Etrans(x::Vector{Float64},g::Grids.DualPatch,ddf_fcn)
    """
    construct body to grid regularization between point `x` and dual
    grid `g`
    """

    # find the support of the discrete delta function    
    rmax = 0.0
    while abs(ddf_fcn(rmax))+abs(ddf_fcn(rmax+1e-5)) > 0.0
        rmax += 0.5
    end

    Extrans = spzeros(length(g.qx))
    Eytrans = spzeros(length(g.qy))

    # rescale the given points into local grid indexing
    xscale = (x[1]-g.xmin[1])/g.Δx
    yscale = (x[2]-g.xmin[2])/g.Δx

    # x edges
    gxscale = [i-0.5 for i=1:size(g.qx,1), j=1:size(g.qx,2)] 
    gyscale = [j for i=1:size(g.qx,1), j=1:size(g.qx,2)] 

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Extrans[ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    # y edges
    gxscale = [i for i=1:size(g.qy,1), j=1:size(g.qy,2)] 
    gyscale = [j-0.5 for i=1:size(g.qy,1), j=1:size(g.qy,2)] 

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eytrans[ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    Extrans,Eytrans


end

function Etrans(b::Bodies.Body,g::Grids.DualPatch,ddf_fcn)
    """
    construct body to grid regularization between body `b` and dual
    grid `g`
    """

    Extrans = spzeros(length(g.qx),b.N)
    Eytrans = spzeros(length(g.qy),b.N)

    for i = 1:b.N
        Etmpx,Etmpy = Etrans(b.x[i],g,ddf_fcn)
    	Extrans[:,i] = Etmpx
	Eytrans[:,i] = Etmpy
    end

    Extrans,Eytrans

end

function Etrans!(dom::Domain)
    dom.Extrans = spzeros(length(dom.grid[1].qx),dom.nbodypts)
    dom.Eytrans = spzeros(length(dom.grid[1].qy),dom.nbodypts)
	 
    for i = 1:dom.nbody
    	Etmpx, Etmpy = Etrans(dom.body[i],dom.grid[1],dom.ddf_fcn)
	dom.Extrans[:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Etmpx
	dom.Eytrans[:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Etmpy
    	 
    end
end

function ECtrans!(dom::Domain)
    dom.ExCtrans = spzeros(length(dom.grid[1].w),dom.nbodypts)
    dom.EyCtrans = spzeros(length(dom.grid[1].w),dom.nbodypts)

    Etrans!(dom)
    qy = zeros(dom.grid[1].qy)
    for i = 1:dom.nbodypts
        qx = reshape(dom.Extrans[:,i],size(dom.grid[1].qx))
    	wtmp = sparse(reshape(Grids.curl(dom.grid[1],qx,qy),length(dom.grid[1].w),1))
	dom.ExCtrans[:,i] = wtmp
    end
    qx = zeros(dom.grid[1].qx)
    for i = 1:dom.nbodypts
        qy = reshape(dom.Eytrans[:,i],size(dom.grid[1].qy))
    	wtmp = sparse(reshape(Grids.curl(dom.grid[1],qx,qy),length(dom.grid[1].w),1))
	dom.EyCtrans[:,i] = wtmp
    end

end


function Base.show(io::IO, dom::Domain)
    println(io, "Domain: xmin = $(dom.xmin), xmax = $(dom.xmax), "*
    		"number of grids = $(dom.ngrid), "*
    		"number of bodies = $(dom.nbody)") 
end



end
