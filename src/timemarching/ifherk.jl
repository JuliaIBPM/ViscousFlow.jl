# IFHERK

"""
    IFHERK(u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[issymmetric=false],[rk::RKParams=RK31])

Construct an integrator to advance a system of the form

du/dt - Au = -B₁ᵀf + r₁(u,t)
B₂u = r₂(u,t)

The resulting integrator will advance the system `(u,f)` by one time step, `Δt`.

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` (by left multiplication) and return same type as `u`
- `B₁ᵀ` : operator acting on type `f` and returning type `u`
- `B₂` : operator acting on type `u` and returning type `f`
- `r₁` : operator acting on type `u` and `t` and returning `u`
- `r₂` : operator acting on type `u` and `t` and returning type `f`
"""
struct IFHERK{NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF}

  # time step size
  Δt :: Float64

  rk :: RKParams

  # Integrating factors
  H :: Vector{FH}

  B₁ᵀ :: FB1  # operates on TF and returns TU
  B₂ :: FB2   # operates on TU and returns TF
  r₁ :: FR1  # function of u and t, returns TU
  r₂ :: FR2  # function of t, returns TF
  P :: FP  # conditioner acting on TF and returning TF

  # Saddle-point systems
  S :: Vector{FS}  # -B₂HB₁ᵀ

  # scratch space
  qᵢ :: TU
  ubuffer :: TU   # should not need this
  w :: Vector{TU}
  fbuffer :: TF

  # flags
  _issymmetric :: Bool

end


function (::Type{IFHERK})(u::TU,f::TF,Δt::Float64,
                          plan_intfact::FI,
                          sys::Tuple{FB1,FB2},
                          rhs::Tuple{FR1,FR2};
                          issymmetric::Bool=false,
                          conditioner::FP = x -> x,
                          store::Bool=false,
                          rk::RKParams{NS}=RK31) where {TU,TF,FI,FB1,FB2,FR1,FR2,FP,NS}


   # templates for the operators
   # B₁ᵀ acts on type TF
   # B₂ acts on TU
   # r₁ acts on TU and time
   # r₂ acts on TU and time
   optypes = ((TF,),(TU,))
   opnames = ("B₁ᵀ","B₂")
   ops = []

   # check for methods for B₁ᵀ and B₂
   for (i,typ) in enumerate(optypes)
     if TU <: Tuple
       opsi = ()
       for I in eachindex(sys[i])
         typI = (typ[1].parameters[I],)
         if method_exists(sys[i][I],typI)
           opsi = (opsi...,sys[i][I])
         elseif method_exists(*,(typeof(sys[i][I]),typI...))
           # generate a method that acts on TU
           opsi = (opsi...,x->sys[i][I]*x)
         else
           error("No valid operator for $(opnames[i]) supplied")
         end
       end
       push!(ops,opsi)
     else
       if method_exists(sys[i],typ)
         push!(ops,sys[i])
       elseif method_exists(*,(typeof(sys[i]),typ...))
         # generate a method that acts on TU
         push!(ops,x->sys[i]*x)
       else
         error("No valid operator for $(opnames[i]) supplied")
       end
     end
   end
   B₁ᵀ, B₂ = ops

   optypes = ((TU,Float64),(TU,Float64))
   opnames = ("r₁","r₂")
   ops = []
   # check for methods for r₁ and r₂
   for (i,typ) in enumerate(optypes)
     if method_exists(rhs[i],typ)
       push!(ops,rhs[i])
     else
       error("No valid operator for $(opnames[i]) supplied")
     end
   end
   r₁, r₂ = ops


    # scratch space
    qᵢ = deepcopy(u)
    ubuffer = deepcopy(u)
    w = [deepcopy(u) for i = 1:NS] # one extra for last step in tuple form
    fbuffer = deepcopy(f)

    dclist = diff([0;rk.c])

    # construct an array of operators for the integrating factor. Each
    # one can act on data of type `u` and return data of the same type.
    # e.g. we can call Hlist[1]*u to get the result.
    if TU <: Tuple
      (FI <: Tuple && length(plan_intfact) == length(u)) ||
                error("plan_intfact argument must be a tuple")
      Hlist = [map((plan,ui) -> plan(dc*Δt,ui),plan_intfact,u) for dc in unique(dclist)]
      Slist = [map((ui,fi,Hi,B₁ᵀi,B₂i) ->
                  SaddleSystem((ui,fi),(Hi,B₁ᵀi,B₂i),issymmetric=issymmetric,isposdef=true,store=store),
                    u,f,H,B₁ᵀ,B₂) for H in Hlist]
    else
      Hlist = [plan_intfact(dc*Δt,u) for dc in unique(dclist)]
      Slist = [SaddleSystem((u,f),(H,B₁ᵀ,B₂),issymmetric=issymmetric,isposdef=true,store=store) for H in Hlist]
    end


    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]
    S = [Slist[i] for i in indexin(dclist,unique(dclist))]

    htype,_ = typeof(H).parameters
    stype,_ = typeof(S).parameters


    # fuse the time step size into the coefficients for some cost savings
    rkdt = deepcopy(rk)
    rkdt.a .*= Δt
    rkdt.c .*= Δt


    ifherksys = IFHERK{NS,htype,typeof(B₁ᵀ),typeof(B₂),typeof(r₁),typeof(r₂),typeof(conditioner),stype,TU,TF}(Δt,rkdt,
                                H,B₁ᵀ,B₂,r₁,r₂,conditioner,S,
                                qᵢ,ubuffer,w,fbuffer,
                                issymmetric)

    # pre-compile
    ifherksys(0.0,u)

    return ifherksys
end

function Base.show(io::IO, scheme::IFHERK{NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF}) where {NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF}
    println(io, "Order-$NS IF-HERK integrator with")
    println(io, "   State of type $TU")
    println(io, "   Force of type $TF")
    println(io, "   Time step size $(scheme.Δt)")
end

# Advance the IFHERK solution by one time step
# This form works when u is a tuple of state vectors
function (scheme::IFHERK{NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF})(t::Float64,u::TU) where
                          {NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU<:Tuple,TF<:Tuple}
  @get scheme (Δt,rk,H,S,r₁,r₂,qᵢ,w,fbuffer,ubuffer)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  i = 1
  for I in eachindex(u)
    ubuffer[I] .= u[I]
    qᵢ[I] .= u[I]
  end
  f = deepcopy(fbuffer)

  if NS > 1
    # first stage, i = 1
    tᵢ₊₁ = t + rk.c[i]

    w[i] = r₁(u,tᵢ₊₁)
    ftmp = r₂(u,tᵢ₊₁) # r₂ # seems like a memory re-allocation...
    for I in eachindex(u)
      w[i][I] .*= rk.a[i,i]  # gᵢ

      ubuffer[I] .+= w[i][I] # r₁ = qᵢ + gᵢ
      fbuffer[I] .= ftmp[I]
      # could solve this system for the whole tuple, too...
      tmp = S[i][I]\(ubuffer[I],fbuffer[I])  # solve saddle point system
      u[I] .= tmp[1]
      f[I] .= tmp[2]

      # diffuse the scratch vectors
      qᵢ[I] .= H[i][I]*qᵢ[I] # qᵢ₊₁ = H(i,i+1)qᵢ
      w[i][I] .= H[i][I]*w[i][I] # H(i,i+1)gᵢ
    end

    # stages 2 through NS-1
    for i = 2:NS-1
      tᵢ₊₁ = t + rk.c[i]

      w[i] = r₁(u,tᵢ₊₁)
      ftmp = r₂(u,tᵢ₊₁) # r₂
      for I in eachindex(u)
        w[i-1][I] .= (w[i-1][I]-S[i-1][I].A⁻¹B₁ᵀf)./(rk.a[i-1,i-1]) # w(i,i-1)

        w[i][I] .*= rk.a[i,i] # gᵢ

        ubuffer[I] .= qᵢ[I] .+ w[i][I] # r₁
        for j = 1:i-1
          ubuffer[I] .+= rk.a[i,j]*w[j][I] # r₁
        end
        fbuffer[I] .= ftmp[I]
        tmp = S[i][I]\(ubuffer[I],fbuffer[I])  # solve saddle point system
        u[I] .= tmp[1]
        f[I] .= tmp[2]

        # diffuse the scratch vectors
        qᵢ[I] .= H[i][I]*qᵢ[I] # qᵢ₊₁ = H(i,i+1)qᵢ
        for j = 1:i
          w[j][I] .= H[i][I]*w[j][I] # for j = i, this sets H(i,i+1)gᵢ
        end
      end

    end
    i = NS
    for I in eachindex(u)
      w[i-1][I] .= (w[i-1][I]-S[i-1][I].A⁻¹B₁ᵀf)./(rk.a[i-1,i-1]) # w(i,i-1)
    end
  end

  # final stage (assembly)
  t = t + rk.c[i]
  w[i] = r₁(u,t)
  ftmp = r₂(u,t) # r₂
  for I in eachindex(u)
    w[i][I] .*= rk.a[i,i]

    ubuffer[I] .= qᵢ[I] .+ w[i][I] # r₁
    for j = 1:i-1
      ubuffer[I] .+= rk.a[i,j]*w[j][I] # r₁
    end
    fbuffer[I] .= ftmp[I]
    tmp = S[i][I]\(ubuffer[I],fbuffer[I])  # solve saddle point system
    u[I] .= tmp[1]
    f[I] .= tmp[2]
    f[I] ./= rk.a[i,i]
  end

  return t, u, f

end

# Advance the IFHERK solution by one time step
function (scheme::IFHERK{NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF})(t::Float64,u::TU) where {NS,FH,FB1,FB2,FR1,FR2,FP,FS,TU,TF}
  @get scheme (Δt,rk,H,S,r₁,r₂,qᵢ,w,fbuffer,ubuffer)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  i = 1
  ubuffer .= u
  qᵢ .= u

  if NS > 1
    # first stage, i = 1
    tᵢ₊₁ = t + rk.c[i]

    w[i] .= rk.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ
    ubuffer .+= w[i] # r₁ = qᵢ + gᵢ
    fbuffer .= r₂(u,tᵢ₊₁) # r₂
    u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system

    # diffuse the scratch vectors
    qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
    w[i] .= H[i]*w[i] # H(i,i+1)gᵢ

    # stages 2 through NS-1
    for i = 2:NS-1
      tᵢ₊₁ = t + rk.c[i]
      w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)./(rk.a[i-1,i-1]) # w(i,i-1)
      w[i] .= rk.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ
      ubuffer .= qᵢ .+ w[i] # r₁
      for j = 1:i-1
        ubuffer .+= rk.a[i,j]*w[j] # r₁
      end
      fbuffer .= r₂(u,tᵢ₊₁) # r₂
      u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system

      #A_ldiv_B!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system

      # diffuse the scratch vectors
      qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
      for j = 1:i
        w[j] .= H[i]*w[j] # for j = i, this sets H(i,i+1)gᵢ
      end

    end
    i = NS
    w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)./(rk.a[i-1,i-1]) # w(i,i-1)
  end

  # final stage (assembly)
  t = t + rk.c[i]
  ubuffer .= qᵢ .+ rk.a[i,i].*r₁(u,t) # r₁
  for j = 1:i-1
    ubuffer .+= rk.a[i,j]*w[j] # r₁
  end
  fbuffer .= r₂(u,t) # r₂
  u, f = S[i]\(ubuffer,fbuffer)  # solve saddle point system
  #A_ldiv_B!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system
  f ./= rk.a[NS,NS]
  return t, u, f

end
