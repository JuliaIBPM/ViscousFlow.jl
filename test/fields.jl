import ViscousFlow: Fields
using FFTW

struct Not{T}
  idx::T
end
import Base: to_indices, uncolon, tail, _maybetail

@inline to_indices(A, inds, I::Tuple{Not, Vararg{Any}}) =
   (setdiff(uncolon(inds, (:, tail(I)...)), I[1].idx), to_indices(A, _maybetail(inds), tail(I))...)

@testset "Grid Routines" begin

  # size
  nx = 12; ny = 12

  # sample point
  i = 5; j = 7
  cellzero = Nodes(Dual,(nx,ny))
  nodezero = Nodes(Primal,cellzero)
  facezero = Edges(Primal,cellzero)
  dualfacezero = Edges(Dual,cellzero)

  @test typeof(cellzero) <: ScalarGridData

  @test typeof(nodezero) <: ScalarGridData

  @test typeof(facezero) <: VectorGridData

  @test typeof(dualfacezero) <: VectorGridData

  @test typeof(facezero.u) <: XEdges{Primal}
  @test typeof(facezero.v) <: YEdges{Primal}

  @test typeof(dualfacezero.u) <: XEdges{Dual}
  @test typeof(dualfacezero.v) <: YEdges{Dual}


  cellunit = deepcopy(cellzero)
  cellunit[i,j] = 1.0

  nodeunit = deepcopy(nodezero)
  nodeunit[i,j] = 1.0

  facexunit = deepcopy(facezero)
  facexunit.u[i,j] = 1.0

  faceyunit = deepcopy(facezero)
  faceyunit.v[i,j] = 1.0

  dualfacexunit = deepcopy(dualfacezero)
  dualfacexunit.u[i,j] = 1.0

  dualfaceyunit = deepcopy(dualfacezero)
  dualfaceyunit.v[i,j] = 1.0

  @testset "Basic array operations" begin
    w = zero(cellunit)
    w .= cellunit
    @test w[i,j] == 1.0
    q = similar(facexunit)
    q .= facexunit
    @test q.u[i,j] == 1.0
    @test iszero(q.v)
  end

  @testset "Inner products and norms" begin
    w = zero(cellunit)
    i0, j0 = rand(2:nx-1), rand(2:ny-1)
    w[i0,j0] = 1.0
    @test Fields.norm(w)*sqrt((nx-2)*(ny-2)) == 1.0
    w .= 1.0
    @test Fields.norm(w) == 1.0

    @test 2*w == w*2

    p = Nodes(Primal,w)
    p .= 1.0
    @test Fields.norm(p) == 1.0
    p2 = deepcopy(p)
    @test Fields.dot(p,p2) == 1.0
    @test Fields.norm(p-p2) == 0.0

    q = Edges(Dual,w)
    q.u .= 1.0
    q2 = deepcopy(q)
    @test Fields.dot(q,q2) == 1.0

    q = Edges(Primal,w)
    q.u .= 1.0
    q2 = deepcopy(q)
    @test Fields.dot(q,q2) == 1.0

    @test Fields.integrate(w) == 1.0

    @test Fields.integrate(p) == 1.0

    q .= 1
    @test Fields.norm(2*q) == sqrt(8)

  end

  @testset "Dual cell center data Laplacian" begin
    lapcell = laplacian(cellunit)
    @test lapcell[i,j] == -4.0
    lapcell[i,j] = 0.0
    @test lapcell[i+1,j] == lapcell[i-1,j] == lapcell[i,j-1] == lapcell[i,j+1] == 1.0
    lapcell[i+1,j] = lapcell[i-1,j] = lapcell[i,j-1] = lapcell[i,j+1] = 0.0
    @test iszero(lapcell)
  end

  @testset "Dual cell center data curl" begin
    q = curl(cellunit)
    @test q.u[i,j-1] == 1.0 && q.u[i,j] == -1.0
    q.u[i,j-1] = q.u[i,j] = 0.0
    @test iszero(q.u)
    @test q.v[i-1,j] == -1.0 && q.v[i,j] == 1.0
    q.v[i-1,j] = q.v[i,j] = 0.0
    @test iszero(q.v)
  end

  @testset "Dual cell node gradient" begin
    q = grad(nodeunit)
    @test q.u[i,j] == 1.0 && q.u[i+1,j] == -1.0
    q.u[i,j] = q.u[i+1,j] = 0.0
    @test iszero(q.u)
    @test q.v[i,j] == 1.0 && q.v[i,j+1] == -1.0
    q.v[i,j] = q.v[i,j+1] = 0.0
    @test iszero(q.v)
  end

  @testset "Face data curl" begin
    cellcurl = curl(facexunit)
    @test cellcurl[i,j] == -1.0 && cellcurl[i,j+1] == 1.0
    cellcurl[i,j] = cellcurl[i,j+1] = 0.0
    @test iszero(cellcurl)
    cellcurl = curl(faceyunit)
    @test cellcurl[i,j] == 1.0 && cellcurl[i+1,j] == -1.0
    cellcurl[i,j] = cellcurl[i+1,j] = 0.0
    @test iszero(cellcurl)
  end

  @testset "Face data divergence" begin
    nodediv = divergence(facexunit)
    @test nodediv[i,j] == -1.0 && nodediv[i-1,j] == 1.0
    nodediv[i,j] = nodediv[i-1,j] = 0.0
    @test iszero(nodediv)
    nodediv = divergence(faceyunit)
    @test nodediv[i,j] == -1.0 && nodediv[i,j-1] == 1.0
    nodediv[i,j] = nodediv[i,j-1] = 0.0
    @test iszero(nodediv)
  end

  @testset "Face data Laplacian" begin
    lap = laplacian(facexunit)
    @test lap.u[i,j] == -4.0
    lap.u[i,j] = 0.0
    @test lap.u[i+1,j] == lap.u[i-1,j] == lap.u[i,j-1] == lap.u[i,j+1] == 1.0
    lap.u[i+1,j] = lap.u[i-1,j] = lap.u[i,j-1] = lap.u[i,j+1] = 0.0
    @test iszero(lap.u)
    @test iszero(lap.v)

    lap = laplacian(faceyunit)
    @test lap.v[i,j] == -4.0
    lap.v[i,j] = 0.0
    @test lap.v[i+1,j] == lap.v[i-1,j] == lap.v[i,j-1] == lap.v[i,j+1] == 1.0
    lap.v[i+1,j] = lap.v[i-1,j] = lap.v[i,j-1] = lap.v[i,j+1] = 0.0
    @test iszero(lap.u)
    @test iszero(lap.v)
  end

  @testset "Dual face data divergence" begin
    celldiv = divergence(dualfacexunit)
    @test celldiv[i,j] == 1.0 && celldiv[i+1,j] == -1.0
    celldiv[i,j] = celldiv[i+1,j] = 0.0
    @test iszero(celldiv)
    celldiv = divergence(dualfaceyunit)
    @test celldiv[i,j] == 1.0 && celldiv[i,j+1] == -1.0
    celldiv[i,j] = celldiv[i,j+1] = 0.0
    @test iszero(celldiv)
  end

  @testset "Face data shift to dual face" begin
    shiftx = Edges(Dual,cellzero)
    cellshift!(shiftx,facexunit)
    @test shiftx.u[i,j] == shiftx.u[i-1,j] == shiftx.u[i,j+1] == shiftx.u[i-1,j+1] == 0.25
    shiftx.u[i,j] = shiftx.u[i-1,j] = shiftx.u[i,j+1] = shiftx.u[i-1,j+1] = 0.0
    @test iszero(shiftx.u)
    @test iszero(shiftx.v)
    shifty = Edges(Dual,cellzero)
    cellshift!(shifty,faceyunit)
    @test shifty.v[i,j] == shifty.v[i,j-1] == shifty.v[i+1,j] == shifty.v[i+1,j-1] == 0.25
    shifty.v[i,j] = shifty.v[i,j-1] = shifty.v[i+1,j] = shifty.v[i+1,j-1] = 0.0
    @test iszero(shifty.u)
    @test iszero(shifty.v)
  end

  @testset "Dual cell center data shift to dual face" begin
    w = Edges(Dual,cellzero)
    cellshift!(w,cellunit)
    @test w.u[i,j] == w.u[i-1,j] == 0.5
    w.u[i,j] = w.u[i-1,j] = 0.0
    @test iszero(w.u)
    @test w.v[i,j] == w.v[i,j-1] == 0.5
    w.v[i,j] = w.v[i,j-1] = 0.0
    @test iszero(w.v)
  end

  @testset "Face data shift to dual cell center" begin
    cellx = Nodes(Dual,cellzero)
    celly = Nodes(Dual,cellzero)
    cellshift!((cellx,celly),facexunit)
    @test cellx[i,j] == 0.5 && cellx[i,j+1] == 0.5
    cellx[i,j] = cellx[i,j+1] = 0.0
    @test iszero(cellx)
    @test iszero(celly)

    cellx = Nodes(Dual,cellzero)
    celly = Nodes(Dual,cellzero)
    cellshift!((cellx,celly),faceyunit)
    @test celly[i,j] == 0.5 && celly[i+1,j] == 0.5
    celly[i,j] = celly[i+1,j] = 0.0
    @test iszero(cellx)
    @test iszero(celly)
  end

  @testset "div curl" begin
    @test iszero(divergence(curl(cellunit)))
    @test iszero(curl(grad(nodeunit)))
  end

  L = plan_laplacian(nx,ny;with_inverse=true)

  @testset "Laplacian of the LGF" begin
    ψ = L\cellunit
    lapψ = L*ψ
    @test lapψ[i,j]≈1.0
    @test isapprox(maximum(abs.(lapψ[Not(i),:])),0.0;atol=10.0*eps()) &&
            isapprox(maximum(abs.(lapψ[:,Not(j)])),0.0;atol=10.0*eps())
  end

end

@testset "Fields" begin
    @testset "Hadamard Product" begin
        edges_p  = Edges{Primal, 30, 40}()
        edges_p.u .= rand(Float64,size(edges_p.u))

        # Should be safe for the output to be the same as the input
        edges_p2 = edges_p ∘ edges_p
        product!(edges_p, edges_p, edges_p)
        @test edges_p2.u == edges_p.u
        @test edges_p2.v == edges_p.v

        edges_d  = Edges{Dual, 30, 40}()
        @test_throws MethodError (edges_p ∘ edges_d)
    end

    @testset "Discrete Laplacian" begin
        s = Nodes{Dual, 30, 40}()
        s[3:end-2, 3:end-2] .= rand(26, 36)

        L = plan_laplacian(30, 40)

        @test L*s ≈ -curl(curl(s))

        @test_throws MethodError (L \ s)

        L = plan_laplacian(30, 40, with_inverse = true, fftw_flags = FFTW.PATIENT)
        @test L \ (L*s) ≈ s

        L! = plan_laplacian!(30, 40, with_inverse = true, fftw_flags = FFTW.PATIENT)
        sold = deepcopy(s)
        L! \ (L!*s)
        @test s ≈ sold
    end

    @testset "Integrating factor" begin
        s = Nodes{Dual, 30, 40}()
        s[15,15] = 1.0

        E1 = plan_intfact(1,s)
        E2 = plan_intfact(2,s)

        @test E1*(E1*s) ≈ E2*s

        E! = plan_intfact!(2,s)
        s2 = deepcopy(s)
        E! * s

        @test s ≈ E2*s2

    end

    @testset "Discrete Divergence" begin
        s = Nodes{Dual, 5, 4}()
        s .= rand(5, 4)

        @test iszero(divergence(curl(s)))

        s = Nodes{Primal, 5, 4}()
        q′ = Edges{Primal, 5, 4}()
        q′.u .= reshape(1:15, 5, 3)
        q′.v .= reshape(1:16, 4, 4)

        # Not sure if this is the behavior we want yet
        # Currently, the ghost cells are not affected
        # by the divergence operator
        #s .= 1.0
        divergence!(s, q′)
        @test s == [ 5.0  5.0  5.0
                     5.0  5.0  5.0
                     5.0  5.0  5.0
                     5.0  5.0  5.0 ]
    end

    @testset "Discrete Curl" begin
        s = Nodes{Dual, 5, 4}()
        s .= reshape(1:20, 4, 5)'

        q = curl(s)

        @test q.u == [ 1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0 ]

        @test q.v == [ -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0 ]
    end

    @testset "Shifting Primal Edges to Dual Edges" begin

        q = Edges{Primal, 5, 4}()
        q.u .= reshape(1:15, 5, 3)
        q.v .= reshape(1:16, 4, 4)
        Qq = Edges{Dual, 5, 4}()

        cellshift!(Qq,q)
        @test Qq.u == [ 0.0  4.0  9.0  0.0
                        0.0  5.0  10.0 0.0
                        0.0  6.0  11.0 0.0
                        0.0  7.0  12.0 0.0 ]

        @test Qq.v == [ 0.0  0.0   0.0
                        3.5  7.5  11.5
                        4.5  8.5  12.5
                        5.5  9.5  13.5
                        0.0  0.0   0.0 ]

    end

    @testset "Shifting Dual Edges to Primal Edges" begin

        q = Edges{Dual, 5, 4}()
        q.u .= reshape(1:16, 4, 4)
        q.v .= reshape(1:15, 5, 3)
        v = Edges{Primal, 5, 4}()
        cellshift!(v,q)

        @test v.u == [ 0.0  0.0   0.0
                       3.5  7.5  11.5
                       4.5  8.5  12.5
                       5.5  9.5  13.5
                       0.0  0.0   0.0 ]

        @test v.v == [ 0.0  4.0   9.0  0.0
                        0.0  5.0  10.0  0.0
                        0.0  6.0  11.0  0.0
                        0.0  7.0  12.0  0.0 ]

    end

    @testset "Shifting Dual Nodes to Dual Edges" begin

        w = Nodes{Dual, 5, 4}()
        w .= reshape(1:20, 5, 4)

        Ww = Edges{Dual, 5, 4}()
        cellshift!(Ww,w)

        @test Ww.u == [ 0.0  6.5  11.5  0.0
                        0.0  7.5  12.5  0.0
                        0.0  8.5  13.5  0.0
                        0.0  9.5  14.5  0.0 ]
        @test Ww.v == [ 0.0   0.0   0.0
                        4.5   9.5  14.5
                        5.5  10.5  15.5
                        6.5  11.5  16.5
                        0.0   0.0   0.0 ]
    end

    @testset "Physical grid" begin

        g = PhysicalGrid((-1.0,3.0),(-2.0,3.0),0.02)
        @test size(g) == (202,252)
        @test size(g,1) == 202
        @test size(g,2) == 252
        @test length(g) == 202*252
        @test origin(g) == (51,101)
        @test cellsize(g) == 0.02
        @test limits(g,1) == (-1.0,3.0)
        @test limits(g,2) == (-2.0,3.0)

    end
end
