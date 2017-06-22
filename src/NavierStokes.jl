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

  Whirl2d.ConstrainedSoln(t,w,f)
end

function TwoLevelBodySoln(dom::Systems.DualDomain)
  # Designated for two-level asymptotic solutions

  t = 0.0
  w = [zeros(dom.grid.cell), zeros(dom.grid.cell)]
  f = [zeros(dom.nbodypts,Whirl2d.ndim),zeros(dom.nbodypts,Whirl2d.ndim)]

  Whirl2d.ConstrainedSoln(t,w,f)
end

#=
This is where we define the operators that define the particular system of equations
we wish to solve
=#
function set_operators!(dom,α)
  # Set up the LGF and integrating factor tables
  Grids.lgf_table!(dom.grid);
  Grids.q_table!(dom.grid,α);

  # Compute the body-to-grid operators
  #Systems.construct_schur!(dom)
  Systems.construct_ECᵀ!(dom)

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = reshape(dom.ECᵀ[1]*f[:,1]+dom.ECᵀ[2]*f[:,2],size(dom.grid.cell))
  # B₂ is function that acts upon data of size s.u and returns data of size s.f
  function B₂(u)
    tmp = reshape(Grids.L⁻¹(dom.grid,u),length(u),1)
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
  S⁻¹(f) = dom.S\reshape(f,length(f),1)
  S₀⁻¹(f) = dom.S₀\reshape(f,length(f),1)

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  r₁(s) = zeros(s.u)
  # r₂ is function that acts upon time value and returns data of size s.f
  r₂(t) = zeros(Float64,dom.ngridpts,2)

  return A⁻¹,B₁ᵀ,B₂,S⁻¹,S₀⁻¹,r₁,r₂
end

function herk!(s::Whirl2d.ConstrainedSoln,p::TimeMarching.TimeParams,dom)





end

end
