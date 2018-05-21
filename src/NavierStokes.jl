module NavierStokes

import Whirl
import Whirl:@get
import Whirl.Grids
import Whirl.Bodies
import Whirl.Systems
import Whirl.TimeMarching

# include("Grids.jl")
# using .Grids
#
# include("Bodies.jl")
# using .Bodies
#
# include("Systems.jl")
# using .Systems
#
# include("TimeMarching.jl")
# using .TimeMarching

# include("./Process.jl")
# using .Process
# include("Grids.jl")
# using .Grids

# For problems without a body
function Soln(dom::Systems.DualDomain)

  #t = 0.0

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Auxiliary data. In this case, streamfunction on dual grid
  #ψ = zeros(dom.grid.cell)

  Whirl.Soln(w)
end

# For Navier-Stokes problems with a body
function BodySoln(dom::Systems.DualDomain)

  # State variable. In this case, vorticity on dual grid
  w = zeros(dom.grid.cell)

  # Constraint force. In this case, Lagrange forces at body points
  f = zeros(dom.nbodypts,Whirl.ndim)

  Whirl.ConstrainedSoln(w,f)
end

# Designated for two-level asymptotic solutions
function TwoLevelBodySoln(dom::Systems.DualDomain)

  w = [zeros(dom.grid.cell), zeros(dom.grid.cell)]
  f = [zeros(dom.nbodypts,Whirl.ndim),zeros(dom.nbodypts,Whirl.ndim)]

  Whirl.ConstrainedSoln(w,f)
end

mutable struct PhysParams
  U∞ :: Vector{Float64}
  Re :: Float64
end

# Eventually, this will hold a lot more specification of the body motions
mutable struct StreamingParams

  "Angular frequency"
  Ω :: Float64

  "A/L, ratio of amplitude to length scale"
  ϵ :: Float64

  "Amplitude factor in x direction. Must be between 0 and 1."
  Ax:: Float64

  "Phase in x direction."
  xϕ :: Float64

  "Amplitude factor in y direction. Must be between 0 and 1."
  Ay :: Float64

  "Phase in y direction."
  yϕ :: Float64

  #"Motion function in x and y directions"
  #X :: Function
end

PhysParams() = PhysParams(zeros(Whirl.ndim),0.0)

StreamingParams() = StreamingParams(0.0,0.0,0.0,0.0,0.0,0.0)

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
  kin = Bodies.RigidBodyMotions.Oscillation(p.Ω,p.Ax,p.xϕ,p.Ay,p.yϕ)
  #X(t) = [p.Ax*sin(p.Ω*t+p.xϕ), p.Ay*sin(p.Ω*t+p.yϕ)]
  #U(t) = p.Ω*[p.Ax*cos(p.Ω*t+p.xϕ), p.Ay*cos(p.Ω*t+p.yϕ)]
  Bodies.set_velocity!(b,Bodies.RigidBodyMotions.RigidBodyMotion(kin),1)
  Bodies.set_velocity!(b,Bodies.RigidBodyMotions.RigidBodyMotion(0.0,0.0),2)
end

# Set forms of the convective term
# ∇⋅(ωu)
function N_divwu(g::Grids.DualPatch, L⁻¹::Grids.Convolution, w::Array{Float64,2},U∞::Array{Float64,1})
  @get Grids (curl, shift, dualdiverg)
   vx, vy = shift(g,curl(g,-L⁻¹(w)))
   wx, wy = shift(g,w)
   return dualdiverg(g,(vx+U∞[1]).*wx./g.Δx,(vy+U∞[2]).*wy./g.Δx)
end

# Set functions that apply body-grid operators
CᵀEᵀ(dom::Systems.DualDomain,f::Array{T,2}) where T =
      reshape(dom.CᵀEᵀ[1]*f[:,1]+dom.CᵀEᵀ[2]*f[:,2],size(dom.grid.cell))

function ẼG̃(dom::Systems.DualDomain,qx::Array{T,2},qy::Array{T,2}) where T
  f = [zeros(dom.nbodypts) for i=1:Whirl.ndim, j = 1:Whirl.ndim]
  f[1,1] = dom.G̃ᵀẼᵀ[1,1]'*squeeze(reshape(qx,length(qx),1),2)
  f[1,2] = dom.G̃ᵀẼᵀ[2,1]'*squeeze(reshape(qx,length(qx),1),2)
  f[2,1] = dom.G̃ᵀẼᵀ[1,2]'*squeeze(reshape(qy,length(qy),1),2)
  f[2,2] = dom.G̃ᵀẼᵀ[2,2]'*squeeze(reshape(qy,length(qy),1),2)
  return f
end

function ECL⁻¹!(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,
                w::Array{T,2},ψ::Array{T,2}) where T
  # This version saves the streamfunction it computes
  ψ = -L⁻¹(w)
  tmp = reshape(-ψ,length(ψ),1)
  return [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
end

function ECL⁻¹(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,w::Array{T,2}) where T
  tmp = reshape(L⁻¹(w),length(w),1)
  return [dom.CᵀEᵀ[1]'*tmp dom.CᵀEᵀ[2]'*tmp]
end

function ẼG̃CL⁻¹(dom::Systems.DualDomain,L⁻¹::Grids.Convolution,w::Array{T,2}) where T
  qx,qy = Grids.curl(dom.grid,-L⁻¹(w))
  return ẼG̃(dom,qx,qy)
end

EÊᵀ(dom::Systems.DualDomain,f::Array{T,2}) where T =
      [dom.Eᵀ[1]'*dom.Êᵀ[1]*f[:,1] dom.Eᵀ[2]'*dom.Êᵀ[2]*f[:,2]]

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

    return Schur{Float64,m}(Sfact,result)
end

function (S⁻¹::Schur{T,M})(f::Array{T,2}) where {T,M}

  copy!(S⁻¹.result,f)
  A_ldiv_B!(S⁻¹.result, S⁻¹.Sfact, reshape(f,length(f),1))

  return [S⁻¹.result[1:M] S⁻¹.result[M+1:2*M]]
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

  curl(ψ) = Grids.curl(dom.grid,ψ)

  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s,t) = -N_divwu(dom.grid,L⁻¹,s.u,physparams.U∞)
  #r₁(s,t) = zeros(s.u)

  return Grids.GridOperators(L⁻¹,curl), TimeMarching.Operators(A⁻¹,r₁)
end


#=
set_operators_body sets up the operators for a Navier-Stokes solution with a body
=#
function set_operators_body!(dom,params,tmp...)
  @get Grids (curl, shift, dualdiverg)
  @get dom (grid,nbodypts)

  physparams = params[1]
  α = params[2]

  # Set up the LGF and integrating factor tables
  if length(dom.grid.lgfhat)==0
    println("Setting up LGF table")
    @time Grids.lgf_table!(dom.grid);
  end
  if length(dom.grid.qhat)==0
    println("Setting up integrating factor table")
    @time Grids.q_table!(dom.grid,α);
  end


  # Id is identity
  Id = Grids.Id(dom.grid)

  # A⁻¹ is function that acts upon data of size s.u and returns data of same size
  A⁻¹ = Grids.Q(dom.grid)

  # L⁻¹ is function that acts upon data of size s.u and returns data of same size
  L⁻¹ = Grids.L⁻¹(dom.grid)
  curl(ψ) = Grids.curl(dom.grid,ψ)

  # Compute the body-to-grid operators
  #Systems.construct_schur!(dom)
  if length(dom.CᵀEᵀ[1])==0
    println("Setting up body-to-grid operator")
    @time Systems.construct_CᵀEᵀ!(dom)
  end

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = CᵀEᵀ(dom,f)

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variable in `s`
  B₂(u) = -ECL⁻¹(dom,L⁻¹,u)

  # Smoother of boundary data (optional)
  #P(f) = Systems.BodyId(dom)(f)
  P(f) = EÊᵀ(dom,f)

  # Compute Schur complements and their inverses
  if length(tmp)==0
    println("Computing Schur complements and inverses")
    @time S⁻¹ = Schur(dom.nbodypts,B₁ᵀ,B₂,A⁻¹)
    @time S₀⁻¹ = Schur(dom.nbodypts,B₁ᵀ,B₂,Id)
  else
    S⁻¹ = tmp[1]
    S₀⁻¹ = tmp[2]
  end


  # r₁ is function that acts upon solution structure s and returns data of size s.u
  # In Navier-Stokes problem, this provides the negative of the convective term.
  r₁(s,t) = -N_divwu(dom.grid,L⁻¹,s.u,physparams.U∞)

  # r₂ is function that acts upon time value and returns data of size s.f
  # In Navier-Stokes problem, this specifies the body velocity (minus the free
  # stream).
  # Should allow U∞ to be time-varying function.
  r₂(s,t) = Systems.Ubody(dom,t) - transpose(repmat(physparams.U∞,1,dom.nbodypts))

  return Grids.GridOperators(L⁻¹,curl),
         TimeMarching.ConstrainedOperators(A⁻¹,B₁ᵀ,B₂,P,S⁻¹,S₀⁻¹,r₁,r₂)
end

function set_operators_two_level_body!(dom,params)

  # Set these operators to accept multiple levels of solution
  physparams = params[1]
  α = params[2]
  sparams = params[3]

  gops1, ops1 = set_operators_body!(dom,params)

  @get ops1 (A⁻¹,B₁ᵀ,B₂,P,S⁻¹,S₀⁻¹,r₁,r₂) (A⁻¹1,B₁ᵀ1,B₂1,P1,S⁻¹1,S₀⁻¹1,r₁1,r₂1)
  @get gops1 (L⁻¹,curl) (L⁻¹1,curl1)

  A⁻¹(u) = A⁻¹1.(u)

  L⁻¹(u) = L⁻¹1.(u)

  function curl(ψ)
    ux = []
    uy = []
    for ψi in ψ
      uxi, uyi = curl1(ψi)
      push!(ux,uxi)
      push!(uy,uyi)
    end
    return ux, uy
  end

  # B₁ᵀ is function that acts upon data of size s.f and returns data of size s.u
  B₁ᵀ(f) = B₁ᵀ1.(f)

  # B₂! is function that takes solution array `s` (and acts upon s.u) and returns
  # data of size s.f. It also modifies the auxiliary variables in `s`
  #B₂!(s::Whirl.SolnType) = [-ECL⁻¹!(dom,L⁻¹1,s.u[1],s.ψ[1]), -ECL⁻¹!(dom,L⁻¹1,s.u[2],s.ψ[2])]
  #B₂(u) = [-ECL⁻¹(dom,L⁻¹1,u[1]), -ECL⁻¹(dom,L⁻¹1,u[2])]
  B₂(u) = -ECL⁻¹.(dom,L⁻¹1,u)

  # Smoother of boundary data (optional)
  P(f) = P1.(f)

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
  #   v[2] = -ΔX(t)⋅∇v[1]
  # where ΔX(t) is the unit displacement vector of the body from its
  # mean position
  function r₂(s,t)
    gradv1 = ẼG̃CL⁻¹(dom,L⁻¹1,s.u[1])/dom.grid.Δx
    vel = zeros(Float64,dom.nbodypts,2)
    for i = 1:dom.nbody
        for ir = dom.firstbpt[i]:dom.firstbpt[i]+dom.body[i].N-1
          # This uses x instead of xtilde, since it is meant to be a perturbation
          # about the current position. But this might need to be fixed.
          # Compare with Ubody calculation.
          x, u, a = dom.body[i].motion[1](t)
          ΔX = real(x) #- dom.body[i].x[ir-dom.firstbpt[i]+1][1]
          ΔY = imag(x) #- dom.body[i].x[ir-dom.firstbpt[i]+1][2]
          vel[ir,1] = -gradv1[1,1][ir]*ΔX-gradv1[1,2][ir]*ΔY
          vel[ir,2] = -gradv1[2,1][ir]*ΔX-gradv1[2,2][ir]*ΔY
        end
    end
    [Systems.Ubody(dom,t,1) - transpose(repmat(physparams.U∞,1,dom.nbodypts)),
     vel]
  end

  return Grids.GridOperators(L⁻¹,curl),
         TimeMarching.ConstrainedOperators(A⁻¹,B₁ᵀ,B₂,P,S⁻¹,S₀⁻¹,r₁,r₂)

end

abstract type FieldType end

struct Field{T} <: FieldType

  t :: Float64
  ω :: T
  ψ :: T
  ux :: T
  uy :: T
  uxfull ::T
  uyfull ::T
  # ψcorr :: T

end

struct ConstrainedField{T,K} <: FieldType

  t :: Float64
  ω :: T
  ψ :: T
  ux :: T
  uy :: T
  uxfull ::T
  uyfull ::T
  # ψcorr :: T
  f :: K

end


function evaluateFields(t::Float64,u::Array{Float64,2},g::Grids.DualPatch,
    gops::Grids.GridOperators)#,nT::UnitRange{Int},Δt::Float64)
  # ω, ψ, ux, uy = evaluateFields(u,g)

  @get gops (L⁻¹,curl)

  ψtmp = -L⁻¹(u)
  utmp, vtmp = curl(ψtmp)
  # ψcorr=Grids.computecorr(g,utmp,vtmp,nT,Δt)
  # ψcorr=utmp
  uxfull=utmp
  uyfull=vtmp
  Field(t,u[g.cellint[1],g.cellint[2]]/g.Δx, ψtmp[g.cellint[1],g.cellint[2]]*g.Δx,
            utmp[g.facexint[1],g.facexint[2]], vtmp[g.faceyint[1],g.faceyint[2]],uxfull,uyfull)
  # Field(t,u[g.cellint[1],g.cellint[2]]/g.Δx, ψtmp[g.cellint[1],g.cellint[2]]*g.Δx,
                      # utmp, vtmp)

end

function evaluateFields(t::Float64,u::Vector{Array{Float64,2}},g::Grids.DualPatch,
    gops::Grids.GridOperators)#,nT::UnitRange{Int},Δt::Float64)

  @get gops (L⁻¹,curl)

  ψtmp = -L⁻¹(u)
  utmp, vtmp = curl(ψtmp)

  ω = []
  for ui in u
    push!(ω,ui[g.cellint[1],g.cellint[2]]/g.Δx)
  end

  ψ = []
  for ψi in ψtmp
    push!(ψ,ψi[g.cellint[1],g.cellint[2]]*g.Δx)
  end

  ux = []
  for uxi in utmp
    push!(ux,uxi[g.facexint[1],g.facexint[2]])

  end

  uy = []
  for uyi in vtmp
    push!(uy,uyi[g.faceyint[1],g.faceyint[2]])

  end

  uxfull = []
  for uxfulli in utmp
    push!(uxfull,uxfulli)
  end

  uyfull = []
  for uyfulli in vtmp
    push!(uyfull,uyfulli)
  end

  # uxfull=utmp
  # uyfull=vtmp
  # ψcorr=computecorr(g,utmp,vtmp,nt,Δt)
  # ψcorr=utmp
  # Field(t, ω, ψ, ux, uy)
  Field(t, ω, ψ, ux, uy, uxfull,uyfull)
end


evaluateFields(s::Whirl.SolnType,g::Grids.DualPatch,gops::Grids.GridOperators)=
      evaluateFields(s.t,s.u,g,gops)


function computecorr(g::Grids.DualPatch,ux,uy,nt::UnitRange{Int},Δt::Float64)
      # map ux and uy to the center of the cells
      ucellx1=[]
      ucelly1=[]
          for i in nt
              push!(ucellx1,Grids.dualshiftx(g,ux[i][1]))
              push!(ucelly1,Grids.dualshifty(g,uy[i][1]))
          end
      # compute the first integral
      intucellx1=[]
      intucelly1=[]

          for  i= 1:length(nt)

              ucellxtemp=ucellx1[1:i]

              intucellx=average(ucellxtemp,1:i).*(i-1)*Δt
              # intucellx=trapzavg(ucellxtemp,1:i).*(i-1)*Δt

              ucellytemp=ucelly1[1:i]

              intucelly=average(ucellytemp,1:i).*(i-1)*Δt
              # intucelly=trapzavg(ucellytemp,1:i).*(i-1)*Δt

              push!(intucellx1,intucellx)
              push!(intucelly1,intucelly)

            end

      # compute the cross product
      cross=[]
          for i=1:length(nt)
              crosstemp=(ucellx1[i].*intucelly1[i]-ucelly1[i].*intucellx1[i])
              push!(cross,crosstemp)
          end
      #correction term

      ψcorrtemp= 0.5.*average(cross,1:length(nt))
      # ψcorrtemp= 0.5.*trapzavg(cross,1:length(nt)))

      ψcorr=ψcorrtemp[g.cellint[1],g.cellint[2]]
      return ψcorr
end

function computecorrtrapz(g::Grids.DualPatch,ux,uy,nt::UnitRange{Int},Δt::Float64)
      # map ux and uy to the center of the cells
      ucellx1=[]
      ucelly1=[]
          for i in nt
              push!(ucellx1,Grids.dualshiftx(g,ux[i][1]))
              push!(ucelly1,Grids.dualshifty(g,uy[i][1]))
          end
      # compute the first integral
      intucellx1=[]
      intucelly1=[]

          for  i= 1:length(nt)

              ucellxtemp=ucellx1[1:i]

              # intucellx=average(ucellxtemp,1:i).*(i-1)*Δt
              intucellx=trapzavg(ucellxtemp,1:i).*(i-1)*Δt

              ucellytemp=ucelly1[1:i]

              # intucelly=average(ucellytemp,1:i).*(i-1)*Δt
              intucelly=trapzavg(ucellytemp,1:i).*(i-1)*Δt

              push!(intucellx1,intucellx)
              push!(intucelly1,intucelly)

            end

      # compute the cross product
      cross=[]
          for i=1:length(nt)
              crosstemp=(ucellx1[i].*intucelly1[i]-ucelly1[i].*intucellx1[i])
              push!(cross,crosstemp)
          end
      #correction term

      # ψcorrtemp= 0.5.*average(cross,1:length(nt)
      ψcorrtemp= 0.5.*trapzavg(cross,1:length(nt))

      ψcorr=ψcorrtemp[g.cellint[1],g.cellint[2]]
end

      # function testpsilgfcorr(g::DualPatch,ux,uy,nt::UnitRange{Int},Δt::Float64,ψ)
      # #corrected psi[t][2]
      #     ψlgf=[]
      #     ψcorr=computecorr(g,ux,uy,nt,Δt)[dom.grid.cellint[1],dom.grid.cellint[2]]
      #     for i in nt
      #         ψlgftemp=ψ[i][2]+ψcorr
      #         push!(ψlgf,ψlgftemp)
      #     end
      #    ψlgf
      # end

# evaluateFields(s::Whirl.SolnType,g::Grids.DualPatch,gops::Grids.GridOperators,nT::UnitRange{Int},Δt::Float64) =
#     evaluateFields(s.t,s.u,g,gops,nT,Δt)


# To evaluate the fields at a specific time `teval`
function (h::Vector{T})(teval::Float64) where {T<:Whirl.SolnType}
  dt(y) = abs.(map(x->x.t,y)-teval)
  imin = find(dt(h)) do x x==minimum(dt(h)) end
  if length(imin) > 0
    return h[imin[end]]
  else
    println("Time $(teval) is not present in this history")
  end
end
trapzavg(f::Union{Vector{Vector{T}},Vector{T}},itr::UnitRange) where T =
        mapreduce(x -> x/length(itr), +,0.*f[itr.start], f[itr.start+1:itr.stop-1])+0.5.*(f[itr.start]+ f[itr.stop])/length(itr)


#Old version of average
average(f::Union{Vector{Vector{T}},Vector{T}},itr::UnitRange) where T =
        mapreduce(x -> x/length(itr), +, f[itr])

# average(f::Union{Vector{Vector{T}},Vector{T}},itr::UnitRange) where T =
# if length(itr)==1 # if the length of the integral is 0
# 0.*f[itr.start]
# else
# mapreduce(x -> x/(length(itr)-1), +,0.*f[itr.start], f[itr.start:itr.stop-1])
# #division by N-1 instead of N,
# # and f(N) doesn't appear in the computing of the average for a 0order approximation
# end

force(s::Whirl.ConstrainedSoln{T,K},g::Grids.DualPatch) where {T,K} = sum(s.f,1)*g.Δx^2

function force(h::Vector{Whirl.ConstrainedSoln{T,K}},g::Grids.DualPatch) where {T,K}
    fh = force.(h,g)
    return map(s -> s.t,h), map(f -> f[1], fh), map(f -> f[2], fh)
end

end
