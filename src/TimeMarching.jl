module TimeMarching

import Whirl:@get

export System, Constrained, Unconstrained, IFRK, IFHERK

"Abstract type for a system of ODEs"
abstract type System{C} end

const Constrained = true
const Unconstrained = false

#function A⁻¹ end
#function r₁ end

struct RKParams{N}
  c::Vector{Float64}
  a::Matrix{Float64}
end

using ..SaddlePointSystems

const RK31 = RKParams{3}([0.5, 1.0, 1.0],
                      [1/2        0        0
                       √3/3 (3-√3)/3        0
                       (3+√3)/6    -√3/3 (3+√3)/6])

const Euler = RKParams{1}([1.0],ones(1,1))

# IFRK

"""
    IFRK(u,Δt,plan_intfact,r₁;[rk::RKParams=RK31])

Construct an operator to advance a system of the form

du/dt - Au = r₁(u,t)

The resulting operator will advance the state `u` by one time step, `Δt`.

# Arguments

- `u` : example of state vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` and return same type as `u`
- `r₁` : operator acting on type `u` and `t` and returning `u`
"""
struct IFRK{NS,FH,FR1,TU}

  # time step size
  Δt :: Float64

  rk :: RKParams

  # Integrating factors
  H :: Vector{FH}

  r₁ :: FR1  # function of u and t, returns TU

  # scratch space
  qᵢ :: TU
  w :: Vector{TU}

end


function (::Type{IFRK})(u::TU,Δt::Float64,
                          plan_intfact::FI,r₁::FR1;
                          rk::RKParams{NS}=RK31) where {TU,FI,FR1,NS}

    # scratch space
    qᵢ = deepcopy(u)
    w = [deepcopy(u) for i = 1:NS-1]

    dclist = diff([0;rk.c])

    # construct an array of operators for the integrating factor. Each
    # one can act on data of type `u` and return data of the same type.
    # e.g. we can call Hlist[1](u) to get the result.
    Hlist = [u->plan_intfact(dc*Δt,u)*u for dc in unique(dclist)]

    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]

    htype,_ = typeof(H).parameters

    IFRK{NS,htype,FR1,TU}(Δt,rk,H,r₁,qᵢ,w)
end

function Base.show(io::IO, scheme::IFRK{NS,FH,FR1,TU}) where {NS,FH,FR1,TU}
    println(io, "Order-$NS IF-RK system with")
    println(io, "   State of type $TU")
    println(io, "   Time step size $(scheme.Δt)")
end

# Advance the IFRK solution by one time step
function (scheme::IFRK{NS,FH,FR1,TU})(t::Float64,u::TU) where {NS,FH,FR1,TU}
  @get scheme (Δt,rk,H,r₁,qᵢ,w)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)

  i = 1
  qᵢ .= u

  if NS > 1
    # first stage, i = 1
    tᵢ₊₁ = t + Δt*rk.c[i]

    w[i] .= Δt*rk.a[i,i]*r₁(u,tᵢ₊₁) # gᵢ

    # diffuse the scratch vectors
    qᵢ .= H[i](qᵢ) # qᵢ₊₁ = H(i,i+1)qᵢ
    w[i] .= H[i](w[i]) # H(i,i+1)gᵢ
    u .= qᵢ .+ w[i]

    # stages 2 through NS-1
    for i = 2:NS-1
      tᵢ₊₁ = t + Δt*rk.c[i]
      w[i-1] ./= Δt*rk.a[i-1,i-1] # w(i,i-1)
      w[i] .= Δt*rk.a[i,i]*r₁(u,tᵢ₊₁) # gᵢ

      # diffuse the scratch vectors and assemble
      qᵢ .= H[i](qᵢ) # qᵢ₊₁ = H(i,i+1)qᵢ
      for j = 1:i
         w[j] .= H[i](w[j]) # for j = i, this sets H(i,i+1)gᵢ
      end
      u .= qᵢ .+ w[i] # r₁
      for j = 1:i-1
        u .+= Δt*rk.a[i,j]*w[j]
      end

    end
    i = NS
    w[i-1] ./= Δt*rk.a[i-1,i-1] # w(i,i-1)
  end

  # final stage (assembly)
  t = t + Δt*rk.c[i]
  u .= qᵢ .+ Δt*rk.a[i,i]*r₁(u,t) # r₁
  for j = 1:i-1
    u .+= Δt*rk.a[i,j]*w[j] # r₁
  end
  u .= H[i](u)

  return t, u

end

# IFHERK

"""
    IFHERK(u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[issymmetric=false],[rk::RKParams=RK31])

Construct an operator to advance a system of the form

du/dt - Au = -B₁ᵀf + r₁(u,t)
B₂u = r₂(t)

The resulting operator will advance the system `(u,f)` by one time step, `Δt`.

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` and return same type as `u`
- `B₁ᵀ` : operator acting on type `f` and returning type `u`
- `B₂` : operator acting on type `u` and returning type `f`
- `r₁` : operator acting on type `u` and `t` and returning `u`
- `r₂` : operator acting on `t` and returning type `f`
"""
struct IFHERK{NS,FH,FB1,FB2,FR1,FR2,TU,TF}

  # time step size
  Δt :: Float64

  rk :: RKParams

  # Integrating factors
  H :: Vector{FH}

  B₁ᵀ :: FB1  # operates on TF and returns TU
  B₂ :: FB2   # operates on TU and returns TF
  r₁ :: FR1  # function of u and t, returns TU
  r₂ :: FR2  # function of t, returns TF

  # Saddle-point systems
  S :: Vector{SaddleSystem}  # -B₂HB₁ᵀ

  # scratch space
  qᵢ :: TU
  ubuffer :: TU   # should not need this
  w :: Vector{TU}
  fbuffer :: TF

  # flags
  _issymmetric :: Bool

end


function (::Type{IFHERK})(u::TU,f::TF,Δt::Float64,
                          plan_intfact::FI,B₁ᵀ::FB1,B₂::FB2,r₁::FR1,r₂::FR2;
                          issymmetric::Bool=false,
                          rk::RKParams{NS}=RK31) where {TU,TF,FI,FB1,FB2,FR1,FR2,NS}

    # scratch space
    qᵢ = deepcopy(u)
    ubuffer = deepcopy(u)
    w = [deepcopy(u) for i = 1:NS-1]
    fbuffer = deepcopy(f)

    dclist = diff([0;rk.c])

    # construct an array of operators for the integrating factor. Each
    # one can act on data of type `u` and return data of the same type.
    # e.g. we can call Hlist[1](u) to get the result.
    Hlist = [u -> plan_intfact(dc*Δt,u)*u for dc in unique(dclist)]

    Slist = [SaddleSystem(u,f,H,B₁ᵀ,B₂,issymmetric=issymmetric,isposdef=true) for H in Hlist]

    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]
    S = [Slist[i] for i in indexin(dclist,unique(dclist))]

    htype,_ = typeof(H).parameters

    IFHERK{NS,htype,FB1,FB2,FR1,FR2,TU,TF}(Δt,rk,
                                H,B₁ᵀ,B₂,r₁,r₂,S,
                                qᵢ,ubuffer,w,fbuffer,
                                issymmetric)
end

function Base.show(io::IO, scheme::IFHERK{NS,FH,FB1,FB2,FR1,FR2,TU,TF}) where {NS,FH,FB1,FB2,FR1,FR2,TU,TF}
    println(io, "Order-$NS IF-HERK system with")
    println(io, "   State of type $TU")
    println(io, "   Force of type $TF")
    println(io, "   Time step size $(scheme.Δt)")
end

# Advance the IFHERK solution by one time step
function (scheme::IFHERK{NS,FH,FB1,FB2,FR1,FR2,TU,TF})(t::Float64,u::TU,f::TF) where {NS,FH,FB1,FB2,FR1,FR2,TU,TF}
  @get scheme (Δt,rk,H,S,B₁ᵀ,B₂,r₁,r₂,qᵢ,w,fbuffer,ubuffer)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)

  i = 1
  ubuffer .= u
  qᵢ .= u

  if NS > 1
    # first stage, i = 1
    tᵢ₊₁ = t + Δt*rk.c[i]

    w[i] .= Δt*rk.a[i,i]*r₁(u,tᵢ₊₁) # gᵢ
    ubuffer .+= w[i] # r₁ = qᵢ + gᵢ
    fbuffer .= r₂(tᵢ₊₁) # r₂
    u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system

    # diffuse the scratch vectors
    qᵢ .= H[i](qᵢ) # qᵢ₊₁ = H(i,i+1)qᵢ
    w[i] .= H[i](w[i]) # H(i,i+1)gᵢ

    # stages 2 through NS-1
    for i = 2:NS-1
      tᵢ₊₁ = t + Δt*rk.c[i]
      w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)/(Δt*rk.a[i-1,i-1]) # w(i,i-1)
      w[i] .= Δt*rk.a[i,i]*r₁(u,tᵢ₊₁) # gᵢ
      ubuffer .= qᵢ .+ w[i] # r₁
      for j = 1:i-1
        ubuffer .+= Δt*rk.a[i,j]*w[j] # r₁
      end
      fbuffer .= r₂(tᵢ₊₁) # r₂
      u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system
      #A_ldiv_B!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system

      # diffuse the scratch vectors
      qᵢ .= H[i](qᵢ) # qᵢ₊₁ = H(i,i+1)qᵢ
      for j = 1:i
        w[j] .= H[i](w[j]) # for j = i, this sets H(i,i+1)gᵢ
      end

    end
    i = NS
    w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)/(Δt*rk.a[i-1,i-1]) # w(i,i-1)
  end

  # final stage (assembly)
  t = t + Δt*rk.c[i]
  ubuffer .= qᵢ .+ Δt*rk.a[i,i]*r₁(u,t) # r₁
  for j = 1:i-1
    ubuffer .+= Δt*rk.a[i,j]*w[j] # r₁
  end
  fbuffer .= r₂(t) # r₂
  u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system
  #A_ldiv_B!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system
  f ./= Δt*rk.a[NS,NS]
  return t, u, f

end


end

#=

struct Operators{TA}

  A⁻¹ :: TA
  r₁ :: Function

end

struct ConstrainedOperators{TA,S,S0}

  A⁻¹ :: TA
  B₁ᵀ :: Function
  B₂ :: Function
  P  :: Function  # Smoother of constraint data (or identity)
  S⁻¹ :: S
  S₀⁻¹ :: S0
  r₁ :: Function
  r₂ :: Function

end


struct TimeParams
  Δt::Float64
  rk::RKParams
end

function Base.show(io::IO, p::TimeParams)
    print(io, "Time step size $(p.Δt)")
end
=#
#=
  We seek to advance, from tⁿ to tⁿ+Δt, the solution of equations of the form

  d/dt (Hu) + HB₁ᵀf =  Hr₁(u,t)

  subject to the constraint

  B₂u = r₂(t)

  where H is the integrating factor H(tⁿ+Δt-t). Each stage of the IF-HERK algorithm
  requires the solution of the saddle point problem

   [Aⁱ  B₁ᵀ](u) = (R₁)
   [B₂ 0  ](f̃) = (R₂)

  where Aⁱ = H(t(i)-t(i+1)). Ultimately, this is solved by use of the Schur complement
  Sⁱ = -B₂(Aⁱ)⁻¹B₁ᵀ. Note that only (Aⁱ)⁻¹ is required, not Aⁱ itself.

  We assume that the algorithm here only has two different forms of (Aⁱ)⁻¹ in the
  various stages: H(Δt/2) and H(0) = I (identity). The former is referred to as
  A⁻¹ in the code. The solution at tⁿ is delivered in the structure `s`, and
  returned as an updated version of itself at tⁿ+Δt. The time parameters are
  provided in `p`. The inverses of the Schur complements corresponding to
  H(Δt/2) and to H(0) are provided as operators S⁻¹ and S₀⁻¹, respectively.

  A⁻¹ is function that acts upon data of size s.u and returns data of same size
  B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₂ is a function that acts upon data of size s and returns data of size s.f
  S and S₀ are factorized Schur complement matrices
  r₁ is function that acts upon solution structure s and returns data of size s.u
  r₂ is function that acts upon time value and returns data of size s.f
=#
#=
function ifherk!(s::Whirl.ConstrainedSoln{T,K},p::TimeParams,ops::ConstrainedOperators) where {T,K}
# Advance the solution by one time step
@get p (Δt,rk)
@get ops (A⁻¹,B₁ᵀ,B₂,P,S⁻¹,S₀⁻¹,r₁,r₂)

# first stage
sᵢ = deepcopy(s)
sᵢ₊₁ = deepcopy(sᵢ)
sᵢ₊₁.t = s.t + Δt*rk.c[1]
A⁻¹gᵢ = Δt*rk.a[1][1]*A⁻¹(r₁(sᵢ,sᵢ₊₁.t))
qᵢ₊₁ = A⁻¹(s.u)
sᵢ₊₁.u = qᵢ₊₁ + A⁻¹gᵢ
sᵢ₊₁.f = -P(S⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ,sᵢ₊₁.t)))
A⁻¹B₁ᵀf = A⁻¹(B₁ᵀ(sᵢ₊₁.f))
@. sᵢ₊₁.u -= A⁻¹B₁ᵀf

w = []
for i = 2:rk.nstage-1
  sᵢ = deepcopy(sᵢ₊₁)
  sᵢ₊₁.t = s.t + Δt*rk.c[i]
  push!(w,(A⁻¹gᵢ-A⁻¹B₁ᵀf)/(Δt*rk.a[i-1][i-1]))
  for j = 1:i-1
    w[j] = A⁻¹(w[j])
  end
  A⁻¹gᵢ = Δt*rk.a[i][i]*A⁻¹(r₁(sᵢ,sᵢ₊₁.t))
  qᵢ₊₁ = A⁻¹(qᵢ₊₁)
  @. sᵢ₊₁.u = qᵢ₊₁ + A⁻¹gᵢ
  for j = 1:i-1
    @. sᵢ₊₁.u += Δt*rk.a[i][j]*w[j]
  end
  sᵢ₊₁.f = -P(S⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ,sᵢ₊₁.t)))
  A⁻¹B₁ᵀf = A⁻¹(B₁ᵀ(sᵢ₊₁.f))
  @. sᵢ₊₁.u -= A⁻¹B₁ᵀf
end

# In final stage, A⁻¹ is assumed to be the identity
i = rk.nstage
sᵢ = deepcopy(sᵢ₊₁)
sᵢ₊₁.t = s.t + Δt*rk.c[i]
push!(w,(A⁻¹gᵢ-A⁻¹B₁ᵀf)/(Δt*rk.a[i-1][i-1]))
for j = 1:i-1
  w[j] = w[j]
end
A⁻¹gᵢ = Δt*rk.a[i][i]*r₁(sᵢ,sᵢ₊₁.t)
sᵢ₊₁.u = qᵢ₊₁ + A⁻¹gᵢ
for j = 1:i-1
  @. sᵢ₊₁.u += Δt*rk.a[i][j]*w[j]
end
sᵢ₊₁.f = -P(S₀⁻¹(B₂(sᵢ₊₁.u) - r₂(sᵢ,sᵢ₊₁.t)))
A⁻¹B₁ᵀf = B₁ᵀ(sᵢ₊₁.f)
@. sᵢ₊₁.u -= A⁻¹B₁ᵀf

# Finalize
s = deepcopy(sᵢ₊₁)
s.f /= Δt*rk.a[rk.nstage][rk.nstage]

return s

end

function ifrk!(s₊, s, Δt, rk::RKParams, sys::System{Unconstrained})
    @get sys (A⁻¹g, q, Ñ, w) # scratch space
    resize!(w, rk.nstage-1)

    u₊ = s₊.u
    u₊ .= s.u

    t₊ = s.t + Δt*rk.c[1]

    A⁻¹(A⁻¹g, r₁(Ñ, u₊, t₊, sys), sys)
    A⁻¹g .*= Δt*rk.a[1, 1]

    A⁻¹(q, s.u, sys)

    @. u₊ = q + A⁻¹g

    for i = 2:rk.nstage-1
        t₊ = s.t + Δt*rk.c[i]

        w[i-1] .= A⁻¹g ./ (Δt*rk.a[i-1, i-1])
        for j = 1:i-1
            A⁻¹(w[j], w[j], sys)
        end

        A⁻¹(A⁻¹g, r₁(Ñ, u₊, t₊, sys), sys)
        A⁻¹g .*= Δt*rk.a[i,i]

        A⁻¹(q, q, sys)

        @. u₊ = q + A⁻¹g

        for j = 1:i-1
            @. u₊ += Δt*rk.a[i,j]*w[j]
        end
    end

    # In final stage, A⁻¹ is assumed to be the identity
    i = rk.nstage
    t₊ = s.t + Δt*rk.c[i]
    w[i-1] .= A⁻¹g ./ (Δt*rk.a[i-1,i-1])

    r₁(A⁻¹g, u₊, t₊, sys)
    A⁻¹g .*= Δt*rk.a[i,i]

    @. u₊ = q + A⁻¹g
    for j in 1:i-1
        u₊ .+= Δt*rk.a[i,j].*w[j]
    end

    s₊.t = t₊

    return s₊
end

=#
