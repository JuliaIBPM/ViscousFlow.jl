module NavierStokes

import Whirl2d
import Whirl2d:@get
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.Systems
import Whirl2d.TimeMarching

# For problems without a body
function Soln(dom::Systems.DualDomain)

  #t = 0.0

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Auxiliary data. In this case, streamfunction on dual grid
  #ψ = zeros(dom.grid.cell)

  Whirl2d.Soln(w)
end

# For Navier-Stokes problems with a body
function BodySoln(dom::Systems.DualDomain)

  #t = 0.0

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Constraint force. In this case, Lagrange forces at body points
  f = zeros(dom.nbodypts,Whirl2d.ndim)

  # Auxiliary data. In this case, streamfunction on dual grid
  #ψ = zeros(dom.grid.cell)

  #Whirl2d.ConstrainedSoln{Float64}(t,w,f,ψ)
  Whirl2d.ConstrainedSoln(w,f)
end

# Designated for two-level asymptotic solutions
function TwoLevelBodySoln(dom::Systems.DualDomain)

  #t = 0.0
  w = [zeros(dom.grid.cell), zeros(dom.grid.cell)]
  f = [zeros(dom.nbodypts,Whirl2d.ndim),zeros(dom.nbodypts,Whirl2d.ndim)]
  #ψ = [zeros(dom.grid.cell), zeros(dom.grid.cell)]

  Whirl2d.ConstrainedSoln(w,f)
end

mutable struct PhysParams
  U∞ :: Array{Float64}
  Re :: Float64
end

# Eventually, this will hold a lot more specification of the body motions
mutable struct StreamingParams

  "Angular frequency"
  Ω :: Float64

  "Amplitude in x direction. Must be between 0 and 1."
  Ax:: Float64

  "Phase in x direction."
  xϕ :: Float64

  "Amplitude in y direction. Must be between 0 and 1."
  Ay :: Float64

  "Phase in y direction."
  yϕ :: Float64

  "Motion function in x and y directions"
  X :: Function
end

PhysParams() = PhysParams(zeros(Whirl2d.ndim),0.0)

StreamingParams() = StreamingParams(0.0,0.0,0.0,0.0,0.0,(t)->[0.0,0.0])

StreamingParams(Ω,Ax,xϕ,Ay,yϕ) = StreamingParams(Ω,Ax,xϕ,Ay,yϕ,(t)->[0.0,0.0])

function set_freestream!(p::PhysParams,U∞)
  p.U∞ = U∞
end

function set_freestream(U∞)
  p = PhysParams()
  set_freestream!(p,U∞)
  p
end

function set_Re!(p::PhysParams,Re)
  p.Re = Re
end

function set_Re(Re)
  p = PhysParams()
  set_Re!(p,Re)
  p
end

function set_oscil_motion!(b::Bodies.Body,p::StreamingParams)
  X(t) = [p.Ax*sin(p.Ω*t+p.xϕ), p.Ay*sin(p.Ω*t+p.yϕ)]
  U(t) = p.Ω*[p.Ax*cos(p.Ω*t+p.xϕ), p.Ay*cos(p.Ω*t+p.yϕ)]
  Bodies.set_velocity!(b,(t,xi)->U(t),1)
  Bodies.set_velocity!(b,(t,xi)->[0.0,0.0],2)
  p.X = X
end

# Set forms of the convective term
# ∇⋅(ωu)
# This version computes a fresh streamfunction
function N_divwu(g::Grids.DualPatch, L⁻¹::Grids.Convolution, u::Array{Float64,2},U∞::Array{Float64,1})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,-L⁻¹(u)))
   wx, wy = shift(g,u)
   dualdiverg(g,(vx+U∞[1]).*wx./g.Δx,(vy+U∞[2]).*wy./g.Δx)
end
# This form makes use of the streamfunction already in the solution structure
function N_divwu(g::Grids.DualPatch,s::Whirl2d.SolnType,U∞::Array{Float64,1})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,s.ψ))
   wx, wy = shift(g,s.u)
   dualdiverg(g,(vx+U∞[1]).*wx,(vy+U∞[2]).*wy)/g.Δx
end

# Set functions that apply body-grid operators
CᵀEᵀ(dom::Systems.DualDomain,f::Array{T,2}) where T =
      reshape(dom.CᵀEᵀ[1]*f[:,1]+dom.CᵀEᵀ[2]*f[:,2],size(dom.grid.cell))

function ẼG̃(dom::Systems.DualDomain,qx::Array{T,2},qy::Array{T,2}) where T
  f = [zeros(dom.nbodypts) for i=1:Whirl2d.ndim, j = 1:Whirl2d.ndim]
  f[1,1] = dom.G̃ᵀẼᵀ[1,1]'*squeeze(reshape(qx,length(qx),1),2)
  f[1,2] = dom.G̃ᵀẼᵀ[2,1]'*squeeze(reshape(qx,length(qx),1),2)
  f[2,1] = dom.G̃ᵀẼᵀ[1,2]'*squeeze(reshape(qy,length(qy),1),2)
  f[2,2] = dom.G̃ᵀẼᵀ[2,2]'*squeeze(reshape(qy,length(qy),1),2)
  f
end

function ECL⁻¹!(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,
                u::Array{T,2},ψ::Array{T,2}) where T
  # This version saves the streamfunction it computes
  ψ = -L⁻¹(u)
  tmp = reshape(-ψ,length(ψ),1)
  [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
end

function ECL⁻¹(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,u::Array{T,2}) where T
  tmp = reshape(L⁻¹(u),length(u),1)
  [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
end

function ẼG̃CL⁻¹(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,u::Array{T,2}) where T
  qx,qy = Grids.curl(dom.grid,-L⁻¹(u))
  ẼG̃(dom,qx,qy)
end

struct Schur{T,M}
    Sfact::Factorization{T}
    result::Array{T, 1}
end

function Schur(m,B₁ᵀ,B₂,A⁻¹::Union{Grids.Convolution,Grids.Identity})

    S = zeros(Float64, 2*m, 2*m)
    result = zeros(Float64, 2*m)

    for i = 1:m
      f = zeros(m,2)
      f[i,1] = 1
      S[:,i] = -reshape(B₂(A⁻¹(B₁ᵀ(f))),length(f),1)
    end
    for i = 1:m
      f = zeros(m,2)
      f[i,2] = 1
      S[:,i+m] = -reshape(B₂(A⁻¹(B₁ᵀ(f))),length(f),1)
    end
    Sfact = factorize(S)

    Schur{Float64,m}(Sfact,result)
end

function (S⁻¹::Schur{T,M})(f::Array{T,2}) where {T,M}

  copy!(S⁻¹.result,f)
  A_ldiv_B!(S⁻¹.result, S⁻¹.Sfact, reshape(f,length(f),1))

  [S⁻¹.result[1:M] S⁻¹.result[M+1:2*M]]
end


#=
This is where we define the operators that define the particular system of equations
we wish to solve.
=#
function set_operators!(dom,params)
  @get Grids (curl, shift, dualdiverg)

  physparams = params[1]
  α = params[2]

  # Set up the LGF and integrating factor tables
  Grids.lgf_table!(dom.grid);
  Grids.q_table!(dom.grid,α);

  # A⁻¹ is function that acts upon data of size s.u and returns data of same size
  A⁻¹ = Grids.Q(dom.grid)

  # L⁻¹ is function that acts upon data of size s.u and returns data of same size
  L⁻¹ = Grids.L⁻¹(dom.grid)

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s,t) = -N_divwu(dom.grid,L⁻¹,s.u,physparams.U∞)
  #r₁(s,t) = zeros(s.u)


  return A⁻¹,L⁻¹,r₁
end


#=
set_operators_body sets up the operators for a Navier-Stokes solution with a body
=#
function set_operators_body!(dom,params)
  @get Grids (curl, shift, dualdiverg)
  @get dom (grid,nbodypts)

  physparams = params[1]
  α = params[2]

  # Set up the LGF and integrating factor tables
  Grids.lgf_table!(dom.grid);
  Grids.q_table!(dom.grid,α);

  # Id is identity
  Id = Grids.Id(dom.grid)

  # A⁻¹ is function that acts upon data of size s.u and returns data of same size
  A⁻¹ = Grids.Q(dom.grid)

  # L⁻¹ is function that acts upon data of size s.u and returns data of same size
  L⁻¹ = Grids.L⁻¹(dom.grid)

  # Compute the body-to-grid operators
  #Systems.construct_schur!(dom)
  Systems.construct_CᵀEᵀ!(dom)

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = CᵀEᵀ(dom,f)

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variable in `s`
  B₂!(s::Whirl2d.SolnType) = -ECL⁻¹!(dom,L⁻¹,s.u,s.ψ)
  B₂(u) = -ECL⁻¹(dom,L⁻¹,u)

  # Compute Schur complements and their inverses
  S⁻¹ = Schur(dom.nbodypts,B₁ᵀ,B₂,A⁻¹)
  S₀⁻¹ = Schur(dom.nbodypts,B₁ᵀ,B₂,Id)


  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s,t) = -N_divwu(dom.grid,L⁻¹,s.u,physparams.U∞)

  # r₂ is function that acts upon time value and returns data of size s.f
  # In Navier-Stokes problem, this specifies the body velocity (minus the free
  # stream).
  # Should allow U∞ to be time-varying function.
  r₂(s,t) = Systems.Ubody(dom,t) - transpose(repmat(physparams.U∞,1,dom.nbodypts))

  return A⁻¹,L⁻¹,B₁ᵀ,B₂!,S⁻¹,S₀⁻¹,r₁,r₂
end

function set_operators_two_level_body!(dom,params)

  # Set these operators to accept multiple levels of solution
  physparams = params[1]
  α = params[2]
  sparams = params[3]

  A⁻¹1,L⁻¹1,B₁ᵀ1,B₂1!,S⁻¹1,S₀⁻¹1,r₁1,r₂1 = set_operators_body!(dom,params)

  A⁻¹(u) = A⁻¹1.(u)

  L⁻¹(u) = L⁻¹1.(u)

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = B₁ᵀ1.(f)

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variables in `s`
  #B₂!(s::Whirl2d.SolnType) = [-ECL⁻¹!(dom,L⁻¹1,s.u[1],s.ψ[1]), -ECL⁻¹!(dom,L⁻¹1,s.u[2],s.ψ[2])]
  #B₂(u) = [-ECL⁻¹(dom,L⁻¹1,u[1]), -ECL⁻¹(dom,L⁻¹1,u[2])]
  B₂!(s::Whirl2d.SolnType) = -ECL⁻¹!.(dom,L⁻¹1,s.u,s.ψ)
  B₂(u) = -ECL⁻¹.(dom,L⁻¹1,u)

  # These functions take in data of size s.f and return data of the same size.
  # They produce the inverse of the Schur complement.
  S⁻¹(f::Vector{Array{T,2}}) where T = S⁻¹1.(f)
  S₀⁻¹(f::Vector{Array{T,2}}) where T = S₀⁻¹1.(f)

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In two-level asymptotic problem, there is no convective term at the first
  # asymptotic level, and at the second level it uses the convective term based
  # on the first level solution
  r₁(s,t) = -[zeros(size(s.u[1])), N_divwu(dom.grid,L⁻¹1,s.u[1],physparams.U∞)]

  # r₂ is function that acts upon time value and returns data of size s.f
  # This specifies the body velocity (minus the free stream).
  # At the second asymptotic level, the boundary data are based on
  #   v[2] = -X(t)⋅∇v[1]
  function r₂(s,t)
    gradv1 = ẼG̃CL⁻¹(dom,L⁻¹1,s.u[1])
    vel = zeros(Float64,dom.nbodypts,2)
    for i = 1:dom.nbody
        ir = dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1
        vel[ir,1] = -gradv1[1,1][ir]*sparams[i].X(t)[1]-gradv1[1,2][ir]*sparams[i].X(t)[2]
        vel[ir,2] = -gradv1[2,1][ir]*sparams[i].X(t)[1]-gradv1[2,2][ir]*sparams[i].X(t)[2]
    end
    [Systems.Ubody(dom,t,1) - transpose(repmat(physparams.U∞,1,dom.nbodypts)),
     vel]
  end

  return A⁻¹,L⁻¹,B₁ᵀ,B₂!,S⁻¹,S₀⁻¹,r₁,r₂

end

end
