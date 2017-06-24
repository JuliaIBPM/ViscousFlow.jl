module NavierStokes

import Whirl2d
import Whirl2d:@get
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.Systems
import Whirl2d.TimeMarching

function Soln(dom::Systems.DualDomain)

  t = 0.0

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Auxiliary data. In this case, streamfunction on dual grid
  ψ = zeros(dom.grid.cell)

  Whirl2d.Soln(t,w,ψ)
end

function BodySoln(dom::Systems.DualDomain)

  t = 0.0

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Constraint force. In this case, Lagrange forces at body points
  f = zeros(dom.nbodypts,Whirl2d.ndim)

  # Auxiliary data. In this case, streamfunction on dual grid
  ψ = zeros(dom.grid.cell)

  Whirl2d.ConstrainedSoln(t,w,f,ψ)
end

function TwoLevelBodySoln(dom::Systems.DualDomain)
  # Designated for two-level asymptotic solutions

  t = 0.0
  w = [zeros(dom.grid.cell), zeros(dom.grid.cell)]
  f = [zeros(dom.nbodypts,Whirl2d.ndim),zeros(dom.nbodypts,Whirl2d.ndim)]
  ψ = [zeros(dom.grid.cell), zeros(dom.grid.cell)]

  Whirl2d.ConstrainedSoln(t,w,f,ψ)
end

mutable struct PhysParams
  U∞ :: Array{Float64}
  Re :: Float64
end

PhysParams() = PhysParams(zeros(Whirl2d.ndim),0.0)

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


# Set forms of the convective term
# ∇⋅(ωu)
function N_divwu(g::Grids.DualPatch,u::Array{Float64,2},U∞::Array{Float64,1})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,-Grids.L⁻¹(g,u)))
   wx, wy = shift(g,u)
   dualdiverg(g,(vx+U∞[1]).*wx,(vy+U∞[2]).*wy)/g.Δx
end
function N_divwu(g::Grids.DualPatch,s::Whirl2d.SolnType,U∞::Array{Float64,1})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,-Grids.L⁻¹(g,s.u)))
   #vx, vy = shift(g,curl(g,s.ψ))
   wx, wy = shift(g,s.u)
   dualdiverg(g,(vx+U∞[1]).*wx,(vy+U∞[2]).*wy)/g.Δx
end

CᵀEᵀ(dom::Systems.DualDomain,f) =
      reshape(dom.CᵀEᵀ[1]*f[:,1]+dom.CᵀEᵀ[2]*f[:,2],size(dom.grid.cell))


function ECL⁻¹!(dom,u,ψ)
  # This version remembers the streamfunction as the auxiliary variable in
  # the solution structure
  ψ = -Grids.L⁻¹(dom.grid,u)
  tmp = reshape(-ψ,length(ψ),1)
  [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
end

function ECL⁻¹(dom,u)
  # This version remembers the streamfunction as the auxiliary variable in
  # the solution structure
  tmp = reshape(Grids.L⁻¹(dom.grid,u),length(u),1)
  [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
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
  A⁻¹(u) = Grids.Q(dom.grid,u)

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s) = -N_divwu(dom.grid,s,physparams.U∞)
  #r₁(s) = zeros(s.u)


  return A⁻¹,r₁
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

  # Compute the body-to-grid operators
  #Systems.construct_schur!(dom)
  Systems.construct_CᵀEᵀ!(dom)

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = CᵀEᵀ(dom,f)

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variable in `s`
  B₂!(s::Whirl2d.SolnType) = -ECL⁻¹!(dom,s.u,s.ψ)
  B₂(u) = -ECL⁻¹(dom,u)

  # A⁻¹ is function that acts upon data of size s.u and returns data of same size
  A⁻¹(u) = Grids.Q(dom.grid,u)

  # Compute Schur complements and their inverses
  m = dom.nbodypts
  S = zeros(2*m,2*m)
  S₀ = zeros(2*m,2*m)
  for i = 1:m
    f = zeros(m,2)
    f[i,1] = 1
    S[:,i] = -reshape(B₂(A⁻¹(B₁ᵀ(f))),length(f),1)
    S₀[:,i] = -reshape(B₂(B₁ᵀ(f)),length(f),1)
  end
  for i = 1:m
    f = zeros(m,2)
    f[i,2] = 1
    S[:,i+m] = -reshape(B₂(A⁻¹(B₁ᵀ(f))),length(f),1)
    S₀[:,i+m] = -reshape(B₂(B₁ᵀ(f)),length(f),1)
  end
  dom.S = factorize(S)
  dom.S₀ = factorize(S₀)

  # These functions take in data of size s.f and return data of the same size.
  # They produce the inverse of the Schur complement.
  S⁻¹(f) = reshape(dom.S\reshape(f,length(f),1),size(f))
  S₀⁻¹(f) = reshape(dom.S₀\reshape(f,length(f),1),size(f))

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s) = -N_divwu(dom.grid,s,physparams.U∞)

  # r₂ is function that acts upon time value and returns data of size s.f
  # In Navier-Stokes problem, this specifies the body velocity (minus the free
  # stream).
  # Should allow U∞ to be time-varying function.
  r₂(t) = Systems.Ubody(dom,t) - transpose(repmat(physparams.U∞,1,dom.nbodypts))

  return A⁻¹,B₁ᵀ,B₂!,S⁻¹,S₀⁻¹,r₁,r₂
end

function set_operators_two_level_body!(dom,params)

  # Set these operators to accept the

  A⁻¹1,B₁ᵀ1,B₂1!,S⁻¹1,S₀⁻¹1,r₁1,r₂1 = set_operators_body!(dom,params)

  A⁻¹(u) = [A⁻¹1(u[1]), A⁻¹1(u[2])]

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = [B₁ᵀ1(f[1]), B₁ᵀ1(f[2])]

  B₂!

end

end
