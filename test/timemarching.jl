import Whirl: Fields, TimeMarching
using Fields
using TimeMarching

@testset "Time Marching" begin

  @testset "Basic Unconstrained" begin

  ω = 4
  u₀ = 1.0
  uex(t) = u₀ + sin(ω*t)/ω

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  rk = RK(u,Δt,r₁,rk=TimeMarching.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u = rk(t,u)
  end

  @test norm(uhist-uex.(T)) ≈ 0.0285193767

  end

  @testset "Unconstrained with integrating factor" begin

  α = 0.5
  ω = 4
  u₀ = 1.0
  uex(t) = u₀*exp(-α*t) + (α*(cos(ω*t)-exp(-α*t))+ω*sin(ω*t))/(α^2+ω^2)

  Fields.plan_intfact(t::Float64,u::Vector{Float64}) = exp(-α*t)

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  ifrk = IFRK(u,Δt,plan_intfact,r₁,rk=TimeMarching.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u = ifrk(t,u)
  end
  @test norm(uhist-uex.(T)) ≈ 0.0180654483

  end

  @testset "Basic unconstrained using IFHERK" begin

  ω = 4
  u₀ = 1.0
  uex(t) = u₀ + sin(ω*t)/ω

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  f = Vector{Float64}()
  TimeMarching.r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  TimeMarching.r₂(u::Vector{Float64},t::Float64) = Vector{Float64}()
  TimeMarching.B₁ᵀ(f::Vector{Float64}) = zeros(Float64,1)
  TimeMarching.B₂(u::Vector{Float64}) = Vector{Float64}()
  Fields.plan_intfact(t::Float64,u::Vector{Float64}) = eye(1)
  ifherk = IFHERK(u,f,Δt,plan_intfact,(B₁ᵀ,B₂),(r₁,r₂),rk=TimeMarching.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u,_ = ifherk(t,u)
  end

  @test norm(uhist-uex.(T)) ≈ 0.0285193767

  end

end
