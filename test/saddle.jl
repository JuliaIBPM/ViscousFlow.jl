import Whirl: Fields, SaddlePointSystems
using Fields
using SaddlePointSystems

@testset "Saddle-Point Systems" begin

nx = 130; ny = 130
Lx = 2.0
dx = Lx/(nx-2)
w = Nodes(Dual,(nx,ny))

L = plan_laplacian(size(w),with_inverse=true)
#L⁻¹(w::T) where {T} = L\w

n = 128
θ = linspace(0,2π,n+1)
R = 0.5
xb = 1.0 + R*cos.(θ[1:n])
yb = 1.0 + R*sin.(θ[1:n])
ds = (2π/n)*R
X = VectorData(xb,yb)
f = ScalarData(X)

E = Regularize(X,dx;ddftype=Fields.Roma,issymmetric=true)
Hmat,Emat = RegularizationMatrix(E,f,w)

PS = SaddleSystem((w,f),(x->L\x,Hmat,Emat),issymmetric=true,isposdef=true)

ψb = ScalarData(X)
w = Nodes(Dual,(nx,ny))
ψb .= -(xb-1)
f .= ones(Float64,n)*ds
ψ = Nodes(Dual,w)
ψ,f = PS\(w,ψb)

fex = -2*cos.(θ[1:n])
@test norm(f./ds-fex,Inf) ≈ 0.035288024
@test ψ[nx,65] ≈ -ψ[1,65]
@test ψ[65,ny] ≈ ψ[65,1]

end
