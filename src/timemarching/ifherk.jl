# IFHERK

"""
    IFHERK(u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[tol=1e-3],[issymmetric=false],[rk::RKParams=RK31])

Construct an integrator to advance a system of the form

du/dt - Au = -B₁ᵀf + r₁(u,t)
B₂u = r₂(u,t)

The resulting integrator will advance the system `(u,f)` by one time step, `Δt`.
The optional argument `tol` sets the tolerance of iterative saddle-point solution,
if applicable.

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` (by left multiplication) and return same type as `u`
- `plan_constraints` : constructor to set up the
- `B₁ᵀ` : operator acting on type `f` and returning type `u`
- `B₂` : operator acting on type `u` and returning type `f`
- `r₁` : operator acting on type `u` and `t` and returning `u`
- `r₂` : operator acting on type `u` and `t` and returning type `f`
"""
struct IFHERK{NS,FH,FR1,FR2,FC,FS,TU,TF}

  # time step size
  Δt :: Float64

  rk :: RKParams
  rkdt :: RKParams

  # Integrating factors
  H :: Vector{FH}

  r₁ :: FR1  # function of u and t, returns TU
  r₂ :: FR2  # function of u and t, returns TF

  # function of u and t, returns B₁ᵀ and B₂
  plan_constraints :: FC

  # Vector of saddle-point systems
  S :: Vector{FS}  # -B₂HB₁ᵀ

  # scratch space
  qᵢ :: TU
  ubuffer :: TU
  w :: Vector{TU}
  fbuffer :: TF

  # iterative solution tolerance
  tol :: Float64

  # flags
  _issymmetric :: Bool # is the system matrix symmetric?
  _isstaticconstraints :: Bool  # do the system constraint operators stay unchanged?
  _isstaticmatrix :: Bool # is the upper left-hand matrix static?
  _isstored :: Bool # is the Schur complement matrix (inverse) stored?

end

function (::Type{IFHERK})(u::TU,f::TF,Δt::Float64,
                          plan_intfact::FI,
                          plan_constraints::FC,
                          rhs::Tuple{FR1,FR2};
                          tol::Float64=1e-3,
                          conditioner::FP = x -> x,
                          issymmetric::Bool=false,
                          isstaticconstraints::Bool=true,
                          isstaticmatrix::Bool=true,
                          isstored::Bool=false,
                          rk::RKParams{NS}=RK31) where {TU,TF,FI,FC,FR1,FR2,FP,NS}


   # templates for the operators
   # r₁ acts on TU and time
   # r₂ acts on TU and time
   optypes = ((TU,Float64),(TU,Float64))
   opnames = ("r₁","r₂")
   ops = []
   # check for methods for r₁ and r₂
   for (i,typ) in enumerate(optypes)
     if hasmethod(rhs[i],typ)
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
      # for each unique element of dclist, create an operator for each
      # element of the tuples u and plan_intfact
      Hlist = map(dc -> map((plan,ui) -> plan(dc*Δt,ui),plan_intfact,u),unique(dclist))
      #Hlist = [map((plan,ui) -> plan(dc*Δt,ui),plan_intfact,u) for dc in unique(dclist)]
    else
      Hlist = map(dc -> plan_intfact(dc*Δt,u),unique(dclist))
      #Hlist = [plan_intfact(dc*Δt,u) for dc in unique(dclist)]
    end
    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]

    # preform the saddle-point systems
    # these are overwritten if B₁ᵀ and B₂ vary with time
    Slist = [construct_saddlesys(plan_constraints,Hi,u,f,0.0,
                    tol,issymmetric,isstored,precompile=false)[1] for Hi in Hlist]
    S = [Slist[i] for i in indexin(dclist,unique(dclist))]

    #S,_ = construct_saddlesys(plan_constraints,rk,H,u,f,0.0,tol,issymmetric,isstored,precompile=false)

    htype,_ = typeof(H).parameters
    stype,_ = typeof(S).parameters


    # fuse the time step size into the coefficients for some cost savings
    rkdt = deepcopy(rk)
    rkdt.a .*= Δt
    rkdt.c .*= Δt


    ifherksys = IFHERK{NS,htype,typeof(r₁),typeof(r₂),FC,stype,TU,TF}(Δt,rk,rkdt,
                                H,r₁,r₂,
                                plan_constraints,S,
                                qᵢ,ubuffer,w,fbuffer,
                                tol,issymmetric,isstaticconstraints,isstaticmatrix,isstored)

    # pre-compile
    #ifherksys(0.0,u)

    return ifherksys
end

function Base.show(io::IO, scheme::IFHERK{NS,FH,FR1,FR2,FC,FS,TU,TF}) where {NS,FH,FR1,FR2,FC,FS,TU,TF}
    println(io, "Order-$NS IF-HERK integrator with")
    println(io, "   State of type $TU")
    println(io, "   Force of type $TF")
    println(io, "   Time step size $(scheme.Δt)")
end

# this function will call the plan_constraints function and return the
# saddle point system for a single instance of H, (and B₁ᵀ and B₂)
# plan_constraints should only compute B₁ᵀ and B₂ (and P if needed)
function construct_saddlesys(plan_constraints::FC,H::FH,
                           u::TU,f::TF,t::Float64,tol::Float64,issymmetric::Bool,isstored::Bool;
                           precompile::Bool=false) where {FC,FH,TU,TF}

    sys = plan_constraints(u,t) # sys contains B₁ᵀ and B₂ before fixing them up

    # B₁ᵀ acts on type TF
    # B₂ acts on TU
    optypes = ((TF,),(TU,))
    opnames = ("B₁ᵀ","B₂")
    ops = []
    # check for methods for B₁ᵀ and B₂
    for (i,typ) in enumerate(optypes)
      if TU <: Tuple
        opsi = ()
        for el in eachindex(sys[i])
          typI = (typ[1].parameters[el],)
          if hasmethod(sys[i][el],typI)
            opsi = (opsi...,sys[i][el])
          elseif hasmethod(*,(typeof(sys[i][el]),typI...))
            # generate a method that acts on TU
            opsi = (opsi...,x->sys[i][el]*x)
          else
            error("No valid operator for $(opnames[i]) supplied")
          end
        end
        push!(ops,opsi)
      else
        if hasmethod(sys[i],typ)
          push!(ops,sys[i])
        elseif hasmethod(*,(typeof(sys[i]),typ...))
          # generate a method that acts on TU
          push!(ops,x->sys[i]*x)
        else
          error("No valid operator for $(opnames[i]) supplied")
        end
      end
    end
    B₁ᵀ, B₂ = ops

    if TU <: Tuple
      S = map((ui,fi,Hi,B₁ᵀi,B₂i) ->
                  SaddleSystem((ui,fi),(Hi,B₁ᵀi,B₂i),tol=tol,issymmetric=issymmetric,isposdef=true,store=isstored,precompile=precompile),
                    u,f,H,B₁ᵀ,B₂)
    else
      S = SaddleSystem((u,f),(H,B₁ᵀ,B₂),tol=tol,
                issymmetric=issymmetric,isposdef=true,store=isstored,precompile=precompile)
    end

    return S, ops



end

# Advance the IFHERK solution by one time step
# This form works when u is a tuple of state vectors
function (scheme::IFHERK{NS,FH,FR1,FR2,FC,FS,TU,TF})(t::Float64,u::TU) where
                          {NS,FH,FR1,FR2,FC,FS,TU<:Tuple,TF<:Tuple}
  @get scheme (rk,rkdt,H,plan_constraints,r₁,r₂,qᵢ,w,fbuffer,ubuffer,tol,
                _isstaticconstraints,_issymmetric,_isstored)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  f = deepcopy(fbuffer)

  i = 1
  tᵢ₊₁ = t
  for el in eachindex(u)
    ubuffer[el] .= u[el]
    qᵢ[el] .= u[el]
  end

  if !_isstaticconstraints
    S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
  else
    S = scheme.S[i]
  end

  if NS > 1
    # first stage, i = 1

    w[i] = r₁(u,tᵢ₊₁)
    ftmp = r₂(u,tᵢ₊₁) # r₂ # seems like a memory re-allocation...
    for el in eachindex(u)
      w[i][el] .*= rkdt.a[i,i]  # gᵢ

      ubuffer[el] .+= w[i][el] # r₁ = qᵢ + gᵢ
      fbuffer[el] .= ftmp[el]
      # could solve this system for the whole tuple, too...
      fill!(f[el],0.0)
      ldiv!((u[el],f[el]),S[el],(ubuffer[el],fbuffer[el]))
      ubuffer[el] .= S[el].A⁻¹B₁ᵀf
      #tmp = S[i][el]\(ubuffer[el],fbuffer[el])  # solve saddle point system
      #u[el] .= tmp[1]
      #f[el] .= tmp[2]

      # diffuse the scratch vectors
      qᵢ[el] .= H[i][el]*qᵢ[el] # qᵢ₊₁ = H(i,i+1)qᵢ
      w[i][el] .= H[i][el]*w[i][el] # H(i,i+1)gᵢ
    end
    tᵢ₊₁ = t + rkdt.c[i]



    # stages 2 through NS-1
    for i = 2:NS-1

      if !_isstaticconstraints
        S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
      else
        S = scheme.S[i]
      end

      w[i] = r₁(u,tᵢ₊₁)
      ftmp = r₂(u,tᵢ₊₁) # r₂
      for el in eachindex(u)
        #w[i-1][el] .= (w[i-1][el]-S[i-1][el].A⁻¹B₁ᵀf)./(rkdt.a[i-1,i-1]) # w(i,i-1)
        w[i-1][el] .= (w[i-1][el]-ubuffer[el])./(rkdt.a[i-1,i-1]) # w(i,i-1)

        w[i][el] .*= rkdt.a[i,i] # gᵢ

        ubuffer[el] .= qᵢ[el] .+ w[i][el] # r₁
        for j = 1:i-1
          ubuffer[el] .+= rkdt.a[i,j]*w[j][el] # r₁
        end
        fbuffer[el] .= ftmp[el]
        fill!(f[el],0.0)
        ldiv!((u[el],f[el]),S[el],(ubuffer[el],fbuffer[el]))
        ubuffer[el] .= S[el].A⁻¹B₁ᵀf

        #tmp = S[i][el]\(ubuffer[el],fbuffer[el])  # solve saddle point system
        #u[el] .= tmp[1]
        #f[el] .= tmp[2]

        # diffuse the scratch vectors
        qᵢ[el] .= H[i][el]*qᵢ[el] # qᵢ₊₁ = H(i,i+1)qᵢ
        for j = 1:i
          w[j][el] .= H[i][el]*w[j][el] # for j = i, this sets H(i,i+1)gᵢ
        end
      end
      tᵢ₊₁ = t + rkdt.c[i]

    end
    i = NS
    if !_isstaticconstraints
      S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
    else
      S = scheme.S[i]
    end
    for el in eachindex(u)
      #w[i-1][el] .= (w[i-1][el]-S[i-1][el].A⁻¹B₁ᵀf)./(rkdt.a[i-1,i-1]) # w(i,i-1)
      w[i-1][el] .= (w[i-1][el]-ubuffer[el])./(rkdt.a[i-1,i-1]) # w(i,i-1)

    end
  end

  # final stage (assembly)
  w[i] = r₁(u,tᵢ₊₁)
  ftmp = r₂(u,tᵢ₊₁) # r₂
  for el in eachindex(u)
    w[i][el] .*= rkdt.a[i,i]

    ubuffer[el] .= qᵢ[el] .+ w[i][el] # r₁
    for j = 1:i-1
      ubuffer[el] .+= rkdt.a[i,j]*w[j][el] # r₁
    end
    fbuffer[el] .= ftmp[el]
    fill!(f[el],0.0)
    ldiv!((u[el],f[el]),S[el],(ubuffer[el],fbuffer[el]))
    #tmp = S[i][el]\(ubuffer[el],fbuffer[el])  # solve saddle point system
    #u[el] .= tmp[1]
    #f[el] .= tmp[2]
    f[el] ./= rkdt.a[i,i]
  end
  t = t + rkdt.c[i]


  return t, u, f

end

# Advance the IFHERK solution by one time step
function (scheme::IFHERK{NS,FH,FR1,FR2,FC,FS,TU,TF})(t::Float64,u::TU) where
                      {NS,FH,FR1,FR2,FC,FS,TU,TF}
  @get scheme (rk,rkdt,H,plan_constraints,r₁,r₂,qᵢ,w,fbuffer,ubuffer,tol,
                    _isstaticconstraints,_issymmetric,_isstored)


  f = deepcopy(fbuffer)

  # H[i] corresponds to H(i,i+1) = H((cᵢ - cᵢ₋₁)Δt)
  # Each of the coefficients includes the time step size

  i = 1
  tᵢ₊₁ = t
  ubuffer .= u
  qᵢ .= u

  if !_isstaticconstraints
    S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
  else
    S = scheme.S[i]
  end


  if NS > 1
    # first stage, i = 1

    w[i] .= rkdt.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ
    ubuffer .+= w[i] # r₁ = qᵢ + gᵢ
    fbuffer .= r₂(u,tᵢ₊₁) # r₂
    u, f = S\(ubuffer,fbuffer)  # solve saddle point system
    ubuffer .= S.A⁻¹B₁ᵀf
    tᵢ₊₁ = t + rkdt.c[i]

    # diffuse the scratch vectors
    qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
    w[i] .= H[i]*w[i] # H(i,i+1)gᵢ



    # stages 2 through NS-1
    for i = 2:NS-1
      if !_isstaticconstraints
        S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
      else
        S = scheme.S[i]
      end
      w[i-1] .= (w[i-1]-ubuffer)./(rkdt.a[i-1,i-1]) # w(i,i-1)
      #w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)./(rkdt.a[i-1,i-1]) # w(i,i-1)
      w[i] .= rkdt.a[i,i].*r₁(u,tᵢ₊₁) # gᵢ
      ubuffer .= qᵢ .+ w[i] # r₁
      for j = 1:i-1
        ubuffer .+= rkdt.a[i,j]*w[j] # r₁
      end
      fbuffer .= r₂(u,tᵢ₊₁) # r₂
      u, f = S\(ubuffer,fbuffer)  # solve saddle point system
      ubuffer .= S.A⁻¹B₁ᵀf
      tᵢ₊₁ = t + rkdt.c[i]

      #ldiv!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system

      # diffuse the scratch vectors
      qᵢ .= H[i]*qᵢ # qᵢ₊₁ = H(i,i+1)qᵢ
      for j = 1:i
        w[j] .= H[i]*w[j] # for j = i, this sets H(i,i+1)gᵢ
      end



    end
    i = NS
    if !_isstaticconstraints
      S,_ = construct_saddlesys(plan_constraints,H[i],u,f,tᵢ₊₁,tol,_issymmetric,_isstored)
    else
      S = scheme.S[i]
    end
    #w[i-1] .= (w[i-1]-S[i-1].A⁻¹B₁ᵀf)./(rkdt.a[i-1,i-1]) # w(i,i-1)
    w[i-1] .= (w[i-1]-ubuffer)./(rkdt.a[i-1,i-1]) # w(i,i-1)
  end

  # final stage (assembly)
  ubuffer .= qᵢ .+ rkdt.a[i,i].*r₁(u,tᵢ₊₁) # r₁
  for j = 1:i-1
    ubuffer .+= rkdt.a[i,j]*w[j] # r₁
  end
  fbuffer .= r₂(u,tᵢ₊₁) # r₂
  u, f = S\(ubuffer,fbuffer)  # solve saddle point system
  #ldiv!((u,f),S[i],(ubuffer,fbuffer)) # solve saddle point system
  f ./= rkdt.a[i,i]
  t = t + rkdt.c[i]

  return t, u, f

end
