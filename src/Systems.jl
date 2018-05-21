module Systems

import Whirl
import Whirl:@get
import Whirl.Grids
import Whirl.Bodies
import Whirl.DDF

using ProgressMeter

abstract type Domain end

mutable struct DualDomain <: Domain
    xmin::Vector{Float64}
    xmax::Vector{Float64}

    grid::Grids.DualPatch

    body::AbstractArray{Bodies.Body}
    nbody::Int
    nbodypts::Int
    firstbpt::Vector{Int}

    ddf_fcn::Function

    """
    Eᵀ operator, acting on body Lagrange points and returning data at
    cell faces
    Ẽᵀ operator, acting on body Lagrange points and returning data at cell
    centers and nodes
    Êᵀ operator, a weighted version of Eᵀ operator
    """
    Eᵀ::Array{SparseMatrixCSC{Float64,Int},1}
    Ẽᵀ::Array{SparseMatrixCSC{Float64,Int},1}
    Êᵀ::Array{SparseMatrixCSC{Float64,Int},1}

    """
    (EC)ᵀ operator, acting on body Lagrange points and returning dual grid cell data
    """
    CᵀEᵀ::Array{SparseMatrixCSC{Float64,Int},1}

    """
    (ẼG̃)ᵀ operator, acting on body Lagrange points and returning grid face data
    """
    G̃ᵀẼᵀ::Array{SparseMatrixCSC{Float64,Int},2}

    """
    Schur complements (factorizations)
    S is Schur complement with the integrating factor, S = -ECL⁻¹QCᵀEᵀ
    S₀ is Schur complement without the integrating factor, S₀ = -ECL⁻¹CᵀEᵀ
    """
    #S::Array{Float64,2}
    #S₀::Array{Float64,2}
    #S
    #S₀

end

function DualDomain()
    xmin = zeros(Whirl.ndim)
    xmax = ones(Whirl.ndim)
    grid = Grids.DualPatch(zeros(Int,Whirl.ndim),0.0,xmin)
    nbody = 0
    nbodypts = 0
    firstbpt = []
    Eᵀ = [spzeros(0) for i=1:Whirl.ndim]
    Ẽᵀ = [spzeros(0) for i=1:Whirl.ndim]
    Êᵀ = [spzeros(0) for i=1:Whirl.ndim]
    CᵀEᵀ = [spzeros(0) for i=1:Whirl.ndim]
    G̃ᵀẼᵀ = [spzeros(0) for i=1:Whirl.ndim,j=1:Whirl.ndim]

    #S = []
    #S₀ = []

    ddf_fcn = Whirl.DDF.ddf_roma
    #ddf_fcn = Whirl.DDF.ddf_goza


    DualDomain(xmin,xmax,grid,[],nbody,nbodypts,firstbpt,
           ddf_fcn,Eᵀ,Ẽᵀ,Êᵀ,CᵀEᵀ,G̃ᵀẼᵀ) #,S,S₀)

end

DualDomain(g::Grids.DualPatch) = add_grid(DualDomain(),g)
DualDomain(b::Bodies.Body) = add_body(DualDomain(),b)
DualDomain(xmin,xmax) = set_dims(DualDomain(),xmin,xmax)


function set_dims!(dom::T,
		               xmin::Vector{Float64},xmax::Vector{Float64}) where T<:Domain
    dom.xmin = xmin
    dom.xmax = xmax
end

function set_dims(dom::T,
		              xmin::Vector{Float64},xmax::Vector{Float64}) where T<:Domain
    set_dims!(dom,xmin,xmax)
    dom
end

function add_grid!(dom::T,g::K) where {T<:Domain,K<:Grids.Grid}
    dom.grid = g

    # resize the domain dimensions to allow the grid to fit
    xmin = [min(dom.xmin[j],g.xmin[j]) for j = 1:Whirl.ndim]
    xmax = [max(dom.xmax[j],g.xmax[j]) for j = 1:Whirl.ndim]
    set_dims!(dom,xmin,xmax)

end

function add_grid(dom::T,g::K) where {T<:Domain,K<:Grids.Grid}
    add_grid!(dom,g)
    dom
end

function add_grid(dom::T,Δx::Float64) where T<:Domain
    # create a domain-filling patch
    N = ceil.(Int,(dom.xmax-dom.xmin)/Δx)
    g = Grids.DualPatch(N,Δx,dom.xmin)
    add_grid!(dom,g)
    dom
end

function add_body!(dom::Domain,b::Bodies.Body)
    dom.nbody += 1
    push!(dom.body,b)

    push!(dom.firstbpt,dom.nbodypts+1)
    dom.nbodypts += b.N

    buff = 0.0

    # resize the domain dimensions to allow the body to fit, with a buffer
    bmin, bmax = Bodies.dims(b)
    xmin = [min(dom.xmin[j],bmin[j]-buff) for j = 1:Whirl.ndim]
    xmax = [max(dom.xmax[j],bmax[j]+buff) for j = 1:Whirl.ndim]
    set_dims!(dom,xmin,xmax)
end

function add_body(dom::Domain,b::Bodies.Body)
    add_body!(dom,b)
    dom
end

BodyId(dom::Domain) = Bodies.Identity(zeros(Float64,dom.nbodypts,Whirl.ndim))


# Body-to-grid operations
function construct_Eᵀ(x::Vector{Float64},g::Grids.DualPatch,ddf_fcn)
    """
    construct body to grid regularization between point `x` and dual
    grid `g`
    """

    # find the support of the discrete delta function
    rmax = 0.0
    while abs(ddf_fcn(rmax))+abs(ddf_fcn(rmax+1e-5)) > 0.0
        rmax += 0.5
    end

    Eᵀ = [spzeros(length(g.facex)),spzeros(length(g.facey))]

    # rescale the given points into local grid indexing
    xscale = (x[1]-g.xmin[1])/g.Δx
    yscale = (x[2]-g.xmin[2])/g.Δx

    # x faces
    gxscale = [i+0.5-g.ifirst[1] for i=1:size(g.facex,1), j=1:size(g.facex,2)]
    gyscale = [j+1.0-g.ifirst[2] for i=1:size(g.facex,1), j=1:size(g.facex,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[1][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    # y faces
    gxscale = [i+1.0-g.ifirst[1] for i=1:size(g.facey,1), j=1:size(g.facey,2)]
    gyscale = [j+0.5-g.ifirst[2] for i=1:size(g.facey,1), j=1:size(g.facey,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[2][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    return Eᵀ

end

function construct_Eᵀ_cellnode(x::Vector{Float64},g::Grids.DualPatch,ddf_fcn)
    """
    construct body to grid regularization between point `x` and dual
    grid `g`
    """

    # find the support of the discrete delta function
    rmax = 0.0
    while abs(ddf_fcn(rmax))+abs(ddf_fcn(rmax+1e-5)) > 0.0
        rmax += 0.5
    end

    Eᵀ = [spzeros(length(g.cell)),spzeros(length(g.node))]

    # rescale the given points into local grid indexing
    xscale = (x[1]-g.xmin[1])/g.Δx
    yscale = (x[2]-g.xmin[2])/g.Δx

    # cell centers
    gxscale = [i+0.5-g.ifirst[1] for i=1:size(g.cell,1), j=1:size(g.cell,2)]
    gyscale = [j+0.5-g.ifirst[2] for i=1:size(g.cell,1), j=1:size(g.cell,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[1][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    # cell nodes
    gxscale = [i+1.0-g.ifirst[1] for i=1:size(g.node,1), j=1:size(g.node,2)]
    gyscale = [j+1.0-g.ifirst[2] for i=1:size(g.node,1), j=1:size(g.node,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[2][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    return Eᵀ

end

function construct_Eᵀ(b::Bodies.Body,g::Grids.DualPatch,ddf_fcn)
    """
    construct body to grid regularization between body `b` and dual
    grid `g`
    """

    Eᵀ = [spzeros(length(g.facex),b.N),spzeros(length(g.facey),b.N)]
    Ẽᵀ = [spzeros(length(g.cell),b.N),spzeros(length(g.node),b.N)]
    Êᵀ = deepcopy(Eᵀ)

    @showprogress 1 "Computing Eᵀ..." for i = 1:b.N
        Eᵀtmp = construct_Eᵀ(b.x[i],g,ddf_fcn)
    	  Eᵀ[1][:,i] = Eᵀtmp[1]
	      Eᵀ[2][:,i] = Eᵀtmp[2]
        Eᵀtmp = construct_Eᵀ_cellnode(b.x[i],g,ddf_fcn)
        Ẽᵀ[1][:,i] = Eᵀtmp[1]
        Ẽᵀ[2][:,i] = Eᵀtmp[2]

    end

    Êᵀ[1] = Diagonal(map(x->x>0.0?1.0/x:0.0,sum(Eᵀ[1],2)[:,1]))*Eᵀ[1]
    Êᵀ[2] = Diagonal(map(x->x>0.0?1.0/x:0.0,sum(Eᵀ[2],2)[:,1]))*Eᵀ[2]

    return Eᵀ, Ẽᵀ, Êᵀ

end

function construct_Eᵀ!(dom::DualDomain)
    dom.Eᵀ[1] = spzeros(length(dom.grid.facex),dom.nbodypts)
    dom.Eᵀ[2] = spzeros(length(dom.grid.facey),dom.nbodypts)

    dom.Ẽᵀ[1] = spzeros(length(dom.grid.cell),dom.nbodypts)
    dom.Ẽᵀ[2] = spzeros(length(dom.grid.node),dom.nbodypts)

    dom.Êᵀ = deepcopy(dom.Eᵀ)

    for i = 1:dom.nbody
    	Eᵀ, Ẽᵀ, Êᵀ = construct_Eᵀ(dom.body[i],dom.grid,dom.ddf_fcn)
	    dom.Eᵀ[1][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀ[1]
	    dom.Eᵀ[2][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀ[2]

      dom.Ẽᵀ[1][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Ẽᵀ[1]
	    dom.Ẽᵀ[2][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Ẽᵀ[2]

      dom.Êᵀ[1][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Êᵀ[1]
	    dom.Êᵀ[2][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Êᵀ[2]

    end
end



function construct_CᵀEᵀ!(dom::DualDomain)

    construct_Eᵀ!(dom)
    @get dom (grid,Eᵀ,Ẽᵀ)

    m = dom.nbodypts

    dom.CᵀEᵀ = [spzeros(length(grid.cell),m) for i = 1:Whirl.ndim]


    qy = zeros(grid.facey)
    @showprogress 1 "Computing CᵀEᵀx..." for i = 1:m
        qx = reshape(Eᵀ[1][:,i],size(grid.facex))
	      dom.CᵀEᵀ[1][:,i] = sparse(reshape(Grids.curl(grid,qx,qy),length(grid.cell),1))
    end
    qx = zeros(grid.facex)
    @showprogress 1 "Computing CᵀEᵀy..." for i = 1:m
        qy = reshape(Eᵀ[2][:,i],size(grid.facey))
	      dom.CᵀEᵀ[2][:,i] = sparse(reshape(Grids.curl(grid,qx,qy),length(grid.cell),1))
    end

    # Construct ẼG̃ operator, which maps grid face data to the body points
    dom.G̃ᵀẼᵀ = Array{SparseMatrixCSC{Float64,Int}}(2,2)
    dom.G̃ᵀẼᵀ[1,1] = spzeros(length(grid.facex),m)
    dom.G̃ᵀẼᵀ[2,1] = spzeros(length(grid.facex),m)
    dom.G̃ᵀẼᵀ[1,2] = spzeros(length(grid.facey),m)
    dom.G̃ᵀẼᵀ[2,2] = spzeros(length(grid.facey),m)


    @showprogress 1 "Computing G̃ᵀẼᵀ..." for i = 1:m
        # Lagrange points to nodes
        facex,facey = Grids.grad(grid,reshape(Ẽᵀ[2][:,i],size(grid.node)))
	      dom.G̃ᵀẼᵀ[1,1][:,i] = sparse(reshape(-facex,length(grid.facex),1))
        dom.G̃ᵀẼᵀ[2,2][:,i] = sparse(reshape(-facey,length(grid.facey),1))
        # Lagrange points to cell centers
        facex,facey = Grids.curl(grid,reshape(Ẽᵀ[1][:,i],size(grid.cell)))
        dom.G̃ᵀẼᵀ[1,2][:,i] = sparse(reshape(facey,length(grid.facey),1))
        dom.G̃ᵀẼᵀ[2,1][:,i] = sparse(reshape(-facex,length(grid.facex),1))
    end


end



# Set the current components of the body velocity from all bodies in the system
# Set them at level l
function Ubody(dom::DualDomain,t::Float64,l=1)::Array{Float64,2}
  vel = zeros(Float64,dom.nbodypts,2)
  for i = 1:dom.nbody
      for j = dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1
        x, u, a = dom.body[i].motion[l](t,dom.body[i].xtilde[j-dom.firstbpt[i]+1])
        vel[j,1] = real(u)
        vel[j,2] = imag(u)
      end
  end
  vel
end

function Base.show(io::IO, dom::DualDomain)
    println(io, "Domain: xmin = $(dom.xmin), xmax = $(dom.xmax)\n"*
        "number of bodies = $(dom.nbody)")
    for i = 1:dom.nbody
        println(io,"$(dom.body[i])")
    end
    if dom.grid.N[1] > 0
      println(io,"$(dom.grid)")
    end

end



end
