module Systems

import Whirl2d
import Whirl2d:@get
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.DDF

abstract type Domain end

mutable struct DualDomain <: Domain
    xmin::Vector{Float64}
    xmax::Vector{Float64}

    grid::Grids.DualPatch

    body::AbstractArray{Bodies.Body}
    nbody::Int
    nbodypts::Int
    firstbpt::Vector{Int}

    ddf_fcn

    """
    E^T operator, acting on body Lagrange points and returning data at
    cell faces
    """
    Eᵀ::Array{SparseMatrixCSC{Float64,Int},1}


    """
    (EC)^T operator, acting on body Lagrange points and returning dual grid cell data
    """
    ECᵀ::Array{SparseMatrixCSC{Float64,Int},1}

    """
    Schur complements (factorizations)
    S is Schur complement with the integrating factor, S = -ECL⁻¹QCᵀEᵀ
    S₀ is Schur complement without the integrating factor, S₀ = -ECL⁻¹CᵀEᵀ
    """
    #S::Array{Float64,2}
    #S₀::Array{Float64,2}
    S
    S₀

end

function DualDomain()
    xmin = zeros(Whirl2d.ndim)
    xmax = ones(Whirl2d.ndim)
    grid = Grids.DualPatch(zeros(Int,Whirl2d.ndim),0.0,xmin)
    nbody = 0
    nbodypts = 0
    firstbpt = []
    Eᵀ = [spzeros(0) for i=1:Whirl2d.ndim]
    ECᵀ = [spzeros(0) for i=1:Whirl2d.ndim]

    S = []
    S₀ = []

    ddf_fcn = Whirl2d.DDF.ddf_roma

    DualDomain(xmin,xmax,grid,[],nbody,nbodypts,firstbpt,
           ddf_fcn,Eᵀ,ECᵀ,S,S₀)

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
    xmin = [min(dom.xmin[j],g.xmin[j]) for j = 1:Whirl2d.ndim]
    xmax = [max(dom.xmax[j],g.xmax[j]) for j = 1:Whirl2d.ndim]
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
    xmin = [min(dom.xmin[j],bmin[j]-buff) for j = 1:Whirl2d.ndim]
    xmax = [max(dom.xmax[j],bmax[j]+buff) for j = 1:Whirl2d.ndim]
    set_dims!(dom,xmin,xmax)
end

function add_body(dom::Domain,b::Bodies.Body)
    add_body!(dom,b)
    dom
end

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
    gxscale = [i-0.5 for i=1:size(g.facex,1), j=1:size(g.facex,2)]
    gyscale = [j for i=1:size(g.facex,1), j=1:size(g.facex,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[1][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    # y faces
    gxscale = [i for i=1:size(g.facey,1), j=1:size(g.facey,2)]
    gyscale = [j-0.5 for i=1:size(g.facey,1), j=1:size(g.facey,2)]

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

    for i = 1:b.N
        Eᵀtmp = construct_Eᵀ(b.x[i],g,ddf_fcn)
    	  Eᵀ[1][:,i] = Eᵀtmp[1]
	      Eᵀ[2][:,i] = Eᵀtmp[2]
    end

    return Eᵀ

end

function construct_Eᵀ!(dom::DualDomain)
    dom.Eᵀ[1] = spzeros(length(dom.grid.facex),dom.nbodypts)
    dom.Eᵀ[2] = spzeros(length(dom.grid.facey),dom.nbodypts)

    for i = 1:dom.nbody
    	Eᵀtmp = construct_Eᵀ(dom.body[i],dom.grid,dom.ddf_fcn)
	    dom.Eᵀ[1][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀtmp[1]
	    dom.Eᵀ[2][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀtmp[2]

    end
end

function construct_ECᵀ!(dom::DualDomain)

    construct_Eᵀ!(dom)
    @get dom (grid,Eᵀ)

    m = dom.nbodypts

    dom.ECᵀ = [spzeros(length(dom.grid.cell),m) for i = 1:Whirl2d.ndim]


    qy = zeros(grid.facey)
    for i = 1:m
        qx = reshape(Eᵀ[1][:,i],size(grid.facex))
	      dom.ECᵀ[1][:,i] = sparse(reshape(Grids.curl(grid,qx,qy),length(grid.cell),1))
    end
    qx = zeros(grid.facex)
    for i = 1:m
        qy = reshape(Eᵀ[2][:,i],size(grid.facey))
	      dom.ECᵀ[2][:,i] = sparse(reshape(Grids.curl(grid,qx,qy),length(grid.cell),1))
    end


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
