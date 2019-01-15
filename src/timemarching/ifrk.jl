# IFRK

"""
    IFRK(u,Δt,plan_intfact,r₁;[rk::RKParams=RK31])

Construct an integrator to advance a system of the form

du/dt - Au = r₁(u,t)

The resulting integrator will advance the state `u` by one time step, `Δt`.

# Arguments

- `u` : example of state vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` (by left multiplication) and return same type as `u`
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
                          plan_intfact::FI,rhs::FR1;
                          rk::RKParams{NS}=RK31) where {TU,FI,FR1,NS}

    # check for methods for r₁ and r₂
    if hasmethod(rhs,(TU,Float64))
        r₁ = rhs
    else
        error("No valid operator for r₁ supplied")
    end

    # scratch space
    qᵢ = deepcopy(u)
    w = [deepcopy(u) for i = 1:NS] # one extra for last step in tuple form

    dclist = diff([0;rk.c])

    # construct an array of operators for the integrating factor. Each
    # one can act on data of type `u` and return data of the same type.
    # e.g. we can call Hlist[1]*u to get the result.
    if TU <: Tuple
      (FI <: Tuple && length(plan_intfact) == length(u)) ||
                error("plan_intfact argument must be a tuple")
      Hlist = map(dc -> map((plan,ui) -> plan(dc*Δt,ui),plan_intfact,u),unique(dclist))                
    else
      Hlist = map(dc -> plan_intfact(dc*Δt,u),unique(dclist))
    end

    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]

    htype,_ = typeof(H).parameters

    # fuse the time step size into the coefficients for some cost savings
    rkdt = deepcopy(rk)
    rkdt.a .*= Δt
    rkdt.c .*= Δt

    ifrksys = IFRK{NS,htype,typeof(r₁),TU}(Δt,rkdt,H,r₁,qᵢ,w)

    #pre-compile
    ifrksys(0.0,u)

    return ifrksys
end

function Base.show(io::IO, scheme::IFRK{NS,FH,FR1,TU}) where {NS,FH,FR1,TU}
    println(io, "Order-$NS IF-RK integator with")
    println(io, "   State of type $TU")
    println(io, "   Time step size $(scheme.Δt)")
end

# Advance the IFRK solution by one time step
# This form works when u is a tuple of state vectors
function (scheme::IFRK{NS,FH,FR1,TU})(t::Float64,u::TU) where {NS,FH,FR1,TU <: Tuple}
  @get scheme (Δt,rk,H,r₁,qᵢ,w)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  i = 1
  tᵢ₊₁ = t
  for I in eachindex(u)
    qᵢ[I] .= u[I]
  end

  if NS > 1
    # first stage, i = 1
    w[i] = r₁(u,tᵢ₊₁)
    for I in eachindex(u)
      w[i][I] .*= rk.a[i,i]  # gᵢ

      # diffuse the scratch vectors
      qᵢ[I] .= H[i][I]*qᵢ[I] # qᵢ₊₁ = H(i,i+1)qᵢ
      w[i][I] .= H[i][I]*w[i][I] # H(i,i+1)gᵢ
      u[I] .= qᵢ[I] .+ w[i][I]
    end
    tᵢ₊₁ = t + rk.c[i]

    # stages 2 through NS-1
    for i = 2:NS-1

      w[i] = r₁(u,tᵢ₊₁)
      for I in eachindex(u)
        w[i-1][I] ./= rk.a[i-1,i-1] # w(i,i-1)

        w[i][I] .*= rk.a[i,i] # gᵢ

        # diffuse the scratch vectors and assemble
        qᵢ[I] .= H[i][I]*qᵢ[I] # qᵢ₊₁ = H(i,i+1)qᵢ
        for j = 1:i
          w[j][I] .= H[i][I]*w[j][I] # for j = i, this sets H(i,i+1)gᵢ
        end
        u[I] .= qᵢ[I] .+ w[i][I] # r₁
        for j = 1:i-1
          u[I] .+= rk.a[i,j]*w[j][I]
        end
      end
      tᵢ₊₁ = t + rk.c[i]


    end
    i = NS
    for I in eachindex(u)
      w[i-1][I] ./= rk.a[i-1,i-1] # w(i,i-1)
    end
  end

  # final stage (assembly)
  w[i] = r₁(u,tᵢ₊₁)
  for I in eachindex(u)
    w[i][I] .*= rk.a[i,i]
    u[I] .= qᵢ[I] .+ w[i][I] # r₁
    for j = 1:i-1
      u[I] .+= rk.a[i,j]*w[j][I] # r₁
    end
    u[I] .= H[i][I]*u[I]
  end
  t = t + rk.c[i]

  return t, u

end

# Advance the IFRK solution by one time step
function (scheme::IFRK{NS,FH,FR1,TU})(t::Float64,u::TU) where {NS,FH,FR1,TU}
  @get scheme (Δt,rk,H,r₁,qᵢ,w)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  i = 1
  tᵢ₊₁ = t
  qᵢ .= u

  if NS > 1
    # first stage, i = 1
    w[i] .= rk.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ

    # diffuse the scratch vectors
    qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
    w[i] .= H[i]*w[i] # H(i,i+1)gᵢ
    u .= qᵢ .+ w[i]
    tᵢ₊₁ = t + rk.c[i]

    # stages 2 through NS-1
    for i = 2:NS-1
      w[i-1] ./= rk.a[i-1,i-1] # w(i,i-1)
      w[i] .= rk.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ

      # diffuse the scratch vectors and assemble
      qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
      for j = 1:i
         w[j] .= H[i]*w[j] # for j = i, this sets H(i,i+1)gᵢ
      end
      u .= qᵢ .+ w[i] # r₁
      for j = 1:i-1
        u .+= rk.a[i,j]*w[j]
      end
      tᵢ₊₁ = t + rk.c[i]

    end
    i = NS
    w[i-1] ./= rk.a[i-1,i-1] # w(i,i-1)
  end

  # final stage (assembly)
  u .= qᵢ .+ rk.a[i,i].*r₁(u,tᵢ₊₁) # r₁
  for j = 1:i-1
    u .+= rk.a[i,j]*w[j] # r₁
  end
  u .= H[i]*u
  t = t + rk.c[i]

  return t, u

end
