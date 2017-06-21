module TimeMarching

import Whirl2d
import Whirl2d:@get

struct RKparams
  nstage::Int
  c::Array{Float64}
  a::Array{Array{Float64}}
end


RK31() = RKparams(3,[0.5,1.0,1.0],
                    [[1/2],[√3/3,(3-√3)/3],[(3+√3)/6,-√3/3,(3+√3)/6]])



struct TimeParams
  Δt::Float64
  rk::RKparams
end

#=
  A⁻¹ is function that acts upon data of size s.u and returns data of same size
  B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₂ is function that acts upon data of size s.u and returns data of size s.f
  S and S₀ are factorized Schur complement matrices
  r₁ is function that acts upon solution structure s and returns data of size s.u
  r₂ is function that acts upon time value and returns data of size s.f
=#
function herk!(s::Whirl2d.Soln,p::TimeParams,A⁻¹,B₁ᵀ,B₂,S⁻¹,S₀⁻¹,r₁,r₂)
# Advance the solution by one time step
@get p (Δt,rk)

sᵢ = deepcopy(s)
sᵢ₊₁ = deepcopy(sᵢ)
sᵢ₊₁.t = s.t + Δt*rk.c[1]
A⁻¹gᵢ = Δt*rk.a[1][1]*A⁻¹(r₁(sᵢ))
qᵢ₊₁ = A⁻¹(s.u)
sᵢ₊₁.u = qᵢ₊₁ + Qgᵢ
sᵢ₊₁.f = -S⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ₊₁.t))
A⁻¹B₁ᵀf = A⁻¹(B₁ᵀ(sᵢ₊₁.f))
sᵢ₊₁.u -= QB₁ᵀf

w = []
for i = 2:rk.nstage
  sᵢ = deepcopy(sᵢ₊₁)
  sᵢ₊₁.t = s.t + Δt*rk.c[i]
  push!(w,(A⁻¹gᵢ-A⁻¹B₁ᵀf)/(Δt*rk.a[i-1][i-1]))
  for j = 1:i-1
    w[j] = A⁻¹(w[j])
  end
  A⁻¹gᵢ = Δt*rk.a[i][i]*A⁻¹(r₁(sᵢ))
  qᵢ₊₁ = A⁻¹(qᵢ₊₁)
  sᵢ₊₁.u = qᵢ₊₁ + A⁻¹gᵢ
  for j = 1:i-1
    sᵢ₊₁.u += Δt*rk.a[i][j]*w[j]
  end
  if i < rk.nstage
    sᵢ₊₁.f = -S⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ₊₁.t))
  else
    sᵢ₊₁.f = -S₀⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ₊₁.t))
  end
  A⁻¹B₁ᵀf = A⁻¹(B₁ᵀ(sᵢ₊₁.f))
  sᵢ₊₁.u -= A⁻¹B₁ᵀf
end

s = deepcopy(sᵢ₊₁)
s.f /= Δt*rk.a[rk.nstage,rk.nstage]

end

end
