module NavierStokes

import Whirl2d
import Whirl2d:@get
import Whirl2d.Grids
import Whirl2d.Bodies
import Whirl2d.Systems
import Whirl2d.TimeMarching

function BodySoln(dom::Systems.DualDomain)

  t = 0.0
  w = zeros(dom.grid.cell)
  f = zeros(dom.nbodypts,Whirl2d.ndim)
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

# Set forms of the convective term
# ∇⋅(ωu)
function N_divwu(g::Grids.DualPatch,u::Array{Float64,2})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = -shift(g,curl(g,Grids.L⁻¹(g,u)))
   wx, wy = shift(g,u)
   dualdiverg(g,(vx+1).*wx,vy.*wy)/g.Δx
end
function N_divwu(g::Grids.DualPatch,s::Whirl2d.Soln)
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,s.ψ))
   wx, wy = shift(g,s.u)
   dualdiverg(g,(vx+1).*wx,vy.*wy)/g.Δx
end

#=
This is where we define the operators that define the particular system of equations
we wish to solve.

set_operators_body sets up the operators for a Navier-Stokes solution with a body
=#
function set_operators_body!(dom,α)
  @get Grids (curl, shift, dualdiverg)
  @get dom (grid,nbodypts)

  # Set up the LGF and integrating factor tables
  Grids.lgf_table!(dom.grid);
  Grids.q_table!(dom.grid,α);

  # Compute the body-to-grid operators
  #Systems.construct_schur!(dom)
  Systems.construct_ECᵀ!(dom)

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = reshape(dom.ECᵀ[1]*f[:,1]+dom.ECᵀ[2]*f[:,2],size(dom.grid.cell))

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variable in `s`
  function B₂!(s::Whirl2d.Soln)
    # This version remembers the streamfunction as the auxiliary variable in
    # the solution structure
    s.ψ = -Grids.L⁻¹(grid,s.u)
    tmp = reshape(s.ψ,length(s.ψ),1)
    [dom.ECᵀ[1]'*tmp dom.ECᵀ[2]'*tmp]
  end
  function B₂(u)
    tmp = -reshape(Grids.L⁻¹(dom.grid,u),length(u),1)
    [dom.ECᵀ[1]'*tmp dom.ECᵀ[2]'*tmp]
  end
  # A⁻¹ is function that acts upon data of size s.u and returns data of same size
  A⁻¹(u) = Grids.Q(dom.grid,u)

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


  # S and S₀ are factorized Schur complement matrices
  S⁻¹(f) = reshape(dom.S\reshape(f,length(f),1),size(f))
  S₀⁻¹(f) = reshape(dom.S₀\reshape(f,length(f),1),size(f))

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  r₁(s) = -N_divwu(dom.grid,s)

  # r₂ is function that acts upon time value and returns data of size s.f
  # Th
  function r₂(t)
    vel = zeros(Float64,dom.nbodypts,2)
    for i = 1:dom.nbody
        for j = dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1
          vel[j,1] = dom.body[i].U(t,dom.body[i].x[j]-dom.firstbpt[i]+1)[1]
          vel[j,2] = dom.body[i].U(t,dom.body[i].x[j]-dom.firstbpt[i]+1)[2]
        end
    end
    vel
  end

  return A⁻¹,B₁ᵀ,B₂!,S⁻¹,S₀⁻¹,r₁,r₂
end

end
