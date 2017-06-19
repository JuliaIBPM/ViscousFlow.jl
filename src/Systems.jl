module Systems

import Whirl2d
import Whirl2d:@get
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.DDF

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

function Domain()
    xmin = zeros(Whirl2d.ndim)
    xmax = ones(Whirl2d.ndim)
    ngrid = nbody = 0
    nbodypts = 0
    firstbpt = []
    Eᵀ = [spzeros(0) for i=1:Whirl2d.ndim]
    ECᵀ = [spzeros(0) for i=1:Whirl2d.ndim]

    S = []
    S₀ = []

    ddf_fcn = Whirl2d.DDF.ddf_roma

    Domain(xmin,xmax,[],ngrid,[],nbody,nbodypts,firstbpt,
           ddf_fcn,Eᵀ,ECᵀ,S,S₀)

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
    push!(dom.grid,g)

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

    Eᵀ = [spzeros(length(g.qx)) for i = 1:Whirl2d.ndim]

    # rescale the given points into local grid indexing
    xscale = (x[1]-g.xmin[1])/g.Δx
    yscale = (x[2]-g.xmin[2])/g.Δx

    # x edges
    gxscale = [i-0.5 for i=1:size(g.qx,1), j=1:size(g.qx,2)]
    gyscale = [j for i=1:size(g.qx,1), j=1:size(g.qx,2)]

    tmpx = abs.(gxscale-xscale).<=rmax
    tmpy = abs.(gyscale-yscale).<=rmax
    ind = find(tmpx.&tmpy)
    Eᵀ[1][ind] = ddf_fcn(gxscale[ind]-xscale).*ddf_fcn(gyscale[ind]-yscale)

    # y edges
    gxscale = [i for i=1:size(g.qy,1), j=1:size(g.qy,2)]
    gyscale = [j-0.5 for i=1:size(g.qy,1), j=1:size(g.qy,2)]

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

    Eᵀ = [spzeros(length(g.qx),b.N) for i = 1:Whirl2d.ndim]

    for i = 1:b.N
        Eᵀtmp = construct_Eᵀ(b.x[i],g,ddf_fcn)
    	  Eᵀ[1][:,i] = Eᵀtmp[1]
	      Eᵀ[2][:,i] = Eᵀtmp[2]
    end

    return Eᵀ

end

function construct_Eᵀ!(dom::Domain)
    dom.Eᵀ[1] = spzeros(length(dom.grid[1].qx),dom.nbodypts)
    dom.Eᵀ[2] = spzeros(length(dom.grid[1].qy),dom.nbodypts)

    for i = 1:dom.nbody
    	Eᵀtmp = construct_Eᵀ(dom.body[i],dom.grid[1],dom.ddf_fcn)
	    dom.Eᵀ[1][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀtmp[1]
	    dom.Eᵀ[2][:,dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1] = Eᵀtmp[2]

    end
end

function construct_ECᵀ!(dom::Domain)

    construct_Eᵀ!(dom)
    @get dom (grid,Eᵀ)

    m = dom.nbodypts

    dom.ECᵀ = [spzeros(length(dom.grid[1].w),m) for i = 1:Whirl2d.ndim]


    qy = zeros(grid[1].qy)
    for i = 1:m
        qx = reshape(Eᵀ[1][:,i],size(grid[1].qx))
	      dom.ECᵀ[1][:,i] = sparse(reshape(Grids.curl(grid[1],qx,qy),length(grid[1].w),1))
    end
    qx = zeros(grid[1].qx)
    for i = 1:m
        qy = reshape(Eᵀ[2][:,i],size(grid[1].qy))
	      dom.ECᵀ[2][:,i] = sparse(reshape(Grids.curl(grid[1],qx,qy),length(grid[1].w),1))
    end


end


function construct_schur!(dom)
  # There is no need to call the ECtrans! routine prior to calling this
  # It is assumed, however, that the grid already has the LGF and integrating
  # factor tables set up

  #=
   Can switch L⁻¹ to L⁻¹_slow and Q to Q_slow below in first set of
  parentheses to evaluate using non-FFT convolution.
  =#
  @get Grids (L⁻¹, Q) (L⁻¹, Q)

  m = dom.nbodypts

  S = zeros(2*m,2*m)
  S₀ = zeros(2*m,2*m)

  # Construct the (EC)^T operator
  construct_ECᵀ!(dom)

  @get dom (grid, ECᵀ)

  # Construct the Schur complement with integrating factor
  #    S = -ECL⁻¹QCᵀEᵀ

  for i = 1:m
    stmpx = Q(grid[1],reshape(ECᵀ[1][:,i],size(grid[1].w)))
    stmpy = Q(grid[1],reshape(ECᵀ[2][:,i],size(grid[1].w)))

    stmpx = reshape(L⁻¹(grid[1],stmpx),length(grid[1].w),1)
    stmpy = reshape(L⁻¹(grid[1],stmpy),length(grid[1].w),1)

    S[1:m,i] = -ECᵀ[1]'*stmpx
    S[1:m,i+m] = -ECᵀ[1]'*stmpy
    S[m+1:2*m,i] = -ECᵀ[2]'*stmpx
    S[m+1:2*m,i+m] = -ECᵀ[2]'*stmpy


  end

  # Factorize the Schur complement
  dom.S = factorize(S);

  # Construct the Schur complement without the integration factor
  #    S₀ = -ECL⁻¹CᵀEᵀ
  for i = 1:m
    stmpx = reshape(
              L⁻¹(grid[1],reshape(ECᵀ[1][:,i],size(grid[1].w))),
              length(grid[1].w),1)
    stmpy = reshape(
              L⁻¹(grid[1],reshape(ECᵀ[2][:,i],size(grid[1].w))),
              length(grid[1].w),1)

    S₀[1:m,i] = -ECᵀ[1]'*stmpx
    S₀[1:m,i+m] = -ECᵀ[1]'*stmpy
    S₀[m+1:2*m,i] = -ECᵀ[2]'*stmpx
    S₀[m+1:2*m,i+m] = -ECᵀ[2]'*stmpy

  end

  # Factorize the Schur complement
  dom.S₀ = factorize(S₀);

end

function Base.show(io::IO, dom::Domain)
    println(io, "Domain: xmin = $(dom.xmin), xmax = $(dom.xmax), "*
    		"number of grids = $(dom.ngrid), "*
    		"number of bodies = $(dom.nbody)")
end



end
