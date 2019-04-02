import ViscousFlow: Fields, SaddlePointSystems

using Compat
using Compat.LinearAlgebra
using Compat: range

@testset "Saddle-Point Systems" begin

nx = 130; ny = 130
Lx = 2.0
dx = Lx/(nx-2)
w = Nodes(Dual,(nx,ny))

L = plan_laplacian(size(w),with_inverse=true)
Linv(x) = L\x
#L⁻¹(w::T) where {T} = L\w

n = 128
θ = range(0,stop=2π,length=n+1)
R = 0.5
xb = 1.0 .+ R*cos.(θ[1:n])
yb = 1.0 .+ R*sin.(θ[1:n])
ds = (2π/n)*R
X = VectorData(xb,yb)
f = ScalarData(X)

E = Regularize(X,dx;ddftype=Fields.Roma,issymmetric=true)
Hmat,Emat = RegularizationMatrix(E,f,w)

PS = SaddleSystem((w,f),(Linv,Hmat,Emat),issymmetric=true,isposdef=true)

ψb = ScalarData(X)
w = Nodes(Dual,(nx,ny))
ψb .= -(xb .- 1)
f .= ones(Float64,n)*ds
ψ = Nodes(Dual,w)
ψ,f = PS\(w,ψb)

fex = -2*cos.(θ[1:n])
@test norm(f./ds-fex,Inf) ≈ 0.035288024
@test ψ[nx,65] ≈ -ψ[1,65]
@test ψ[65,ny] ≈ ψ[65,1]

ru = ones(Float64,2)
rf = Vector{Float64}()
A⁻¹(u::Vector{Float64}) = u
B₁ᵀ(f::Vector{Float64}) = zeros(Float64,2)
B₂(u::Vector{Float64}) = Vector{Float64}()
sys = SaddleSystem((ru,rf),(A⁻¹,B₁ᵀ,B₂))

u,fnull = sys\(ru,rf)
@test u ≈ [1.0,1.0]
@test length(fnull) == 0

sysys = (PS,sys)
rhs = ((w,ψb),(ru,rf))
(ψ,f),(u,fnull) = sysys\rhs

@test ψ[nx,65] ≈ -ψ[1,65]
@test ψ[65,ny] ≈ ψ[65,1]
@test u ≈ [1.0,1.0]
@test length(fnull) == 0


end
