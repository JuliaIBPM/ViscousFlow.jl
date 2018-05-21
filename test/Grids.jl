# This stuff is for complementary indexing
immutable Not{T}
  idx::T
end
import Base: to_indices, uncolon, tail, _maybetail

@inline to_indices(A, inds, I::Tuple{Not, Vararg{Any}}) =
   (setdiff(uncolon(inds, (:, tail(I)...)), I[1].idx), to_indices(A, _maybetail(inds), tail(I))...)


@testset "Grid Routines" begin

  xmin = [-1.0,-1.0]
  xmax = [1.0,1.0]
  dom = Whirl.Systems.DualDomain(xmin,xmax)

  Δx = 0.2
  dom = Whirl.Systems.add_grid(dom,Δx)

  i = 5; j = 7

  cellzero = zeros(dom.grid.cell)
  nodezero = zeros(dom.grid.node)
  facexzero = zeros(dom.grid.facex)
  faceyzero = zeros(dom.grid.facey)
  dualfacexzero = zeros(dom.grid.dualfacex)
  dualfaceyzero = zeros(dom.grid.dualfacex)

  cellunit = copy(cellzero)
  cellunit[i,j] = 1.0

  nodeunit = copy(nodezero)
  nodeunit[i,j] = 1.0

  facexunit = copy(facexzero)
  facexunit[i,j] = 1.0

  faceyunit = copy(faceyzero)
  faceyunit[i,j] = 1.0

  dualfacexunit = copy(dualfacexzero)
  dualfacexunit[i,j] = 1.0

  dualfaceyunit = copy(dualfaceyzero)
  dualfaceyunit[i,j] = 1.0


  @testset "Dual cell center data Laplacian" begin
    lapcell = Whirl.Grids.lap(dom.grid,cellunit)
    @test lapcell[i,j] == -4.0
    lapcell[i,j] = 0.0
    @test lapcell[i+1,j] == lapcell[i-1,j] == lapcell[i,j-1] == lapcell[i,j+1] == 1.0
    lapcell[i+1,j] = lapcell[i-1,j] = lapcell[i,j-1] = lapcell[i,j+1] = 0.0
    @test all(lapcell.==0.0)
  end

  @testset "Dual cell center data curl" begin
    qx,qy = Whirl.Grids.curl(dom.grid,cellunit)
    @test qx[i,j-1] == 1.0 && qx[i,j] == -1.0
    qx[i,j-1] = qx[i,j] = 0.0
    @test all(qx.==0.0)
    @test qy[i-1,j] == -1.0 && qy[i,j] == 1.0
    qy[i-1,j] = qy[i,j] = 0.0
    @test all(qy.==0.0)
  end

  @testset "Dual cell node gradient" begin
    qx,qy = Whirl.Grids.grad(dom.grid,nodeunit)
    @test qx[i,j] == 1.0 && qx[i+1,j] == -1.0
    qx[i,j] = qx[i+1,j] = 0.0
    @test all(qx.==0.0)
    @test qy[i,j] == 1.0 && qy[i,j+1] == -1.0
    qy[i,j] = qy[i,j+1] = 0.0
    @test all(qy.==0.0)
  end

  @testset "Face data curl" begin
    cellcurl = Whirl.Grids.curl(dom.grid,facexunit,faceyzero)
    @test cellcurl[i,j] == -1.0 && cellcurl[i,j+1] == 1.0
    cellcurl[i,j] = cellcurl[i,j+1] = 0.0
    @test all(cellcurl.==0.0)
    cellcurl = Whirl.Grids.curl(dom.grid,facexzero,faceyunit)
    @test cellcurl[i,j] == 1.0 && cellcurl[i+1,j] == -1.0
    cellcurl[i,j] = cellcurl[i+1,j] = 0.0
    @test all(cellcurl.==0.0)
  end

  @testset "Face data divergence" begin
    nodediv = Whirl.Grids.diverg(dom.grid,facexunit,faceyzero)
    @test nodediv[i,j] == -1.0 && nodediv[i-1,j] == 1.0
    nodediv[i,j] = nodediv[i-1,j] = 0.0
    @test all(nodediv.==0.0)
    nodediv = Whirl.Grids.diverg(dom.grid,facexzero,faceyunit)
    @test nodediv[i,j] == -1.0 && nodediv[i,j-1] == 1.0
    nodediv[i,j] = nodediv[i,j-1] = 0.0
    @test all(nodediv.==0.0)
  end

  @testset "Face data Laplacian" begin
    lapx,lapy = Whirl.Grids.lap(dom.grid,facexunit,faceyzero)
    @test lapx[i,j] == -4.0
    lapx[i,j] = 0.0
    @test lapx[i+1,j] == lapx[i-1,j] == lapx[i,j-1] == lapx[i,j+1] == 1.0
    lapx[i+1,j] = lapx[i-1,j] = lapx[i,j-1] = lapx[i,j+1] = 0.0
    @test all(lapx.==0.0)
    @test all(lapy.==0.0)

    lapx,lapy = Whirl.Grids.lap(dom.grid,facexzero,faceyunit)
    @test lapy[i,j] == -4.0
    lapy[i,j] = 0.0
    @test lapy[i+1,j] == lapy[i-1,j] == lapy[i,j-1] == lapy[i,j+1] == 1.0
    lapy[i+1,j] = lapy[i-1,j] = lapy[i,j-1] = lapy[i,j+1] = 0.0
    @test all(lapx.==0.0)
    @test all(lapy.==0.0)
  end

  @testset "Dual face data divergence" begin
    celldiv = Whirl.Grids.dualdiverg(dom.grid,dualfacexunit,dualfaceyzero)
    @test celldiv[i,j] == 1.0 && celldiv[i+1,j] == -1.0
    celldiv[i,j] = celldiv[i+1,j] = 0.0
    @test all(celldiv.==0.0)
    celldiv = Whirl.Grids.dualdiverg(dom.grid,dualfacexzero,dualfaceyunit)
    @test celldiv[i,j] == 1.0 && celldiv[i,j+1] == -1.0
    celldiv[i,j] = celldiv[i,j+1] = 0.0
    @test all(celldiv.==0.0)
  end

  @testset "Face data shift to dual face" begin
    shiftx,shifty = Whirl.Grids.shift(dom.grid,facexunit,faceyzero)
    @test shiftx[i,j] == shiftx[i-1,j] == shiftx[i,j+1] == shiftx[i-1,j+1] == 0.25
    shiftx[i,j] = shiftx[i-1,j] = shiftx[i,j+1] = shiftx[i-1,j+1] = 0.0
    @test all(shiftx.==0.0)
    @test all(shifty.==0.0)
    shiftx,shifty = Whirl.Grids.shift(dom.grid,facexzero,faceyunit)
    @test shifty[i,j] == shifty[i,j-1] == shifty[i+1,j] == shifty[i+1,j-1] == 0.25
    shifty[i,j] = shifty[i,j-1] = shifty[i+1,j] = shifty[i+1,j-1] = 0.0
    @test all(shiftx.==0.0)
    @test all(shifty.==0.0)
  end

  @testset "Dual cell center data shift to dual face" begin
    shiftx,shifty = Whirl.Grids.shift(dom.grid,cellunit)
    @test shiftx[i,j] == shiftx[i-1,j] == 0.5
    shiftx[i,j] = shiftx[i-1,j] = 0.0
    @test all(shiftx.==0.0)
    @test shifty[i,j] == shifty[i,j-1] == 0.5
    shifty[i,j] = shifty[i,j-1] = 0.0
    @test all(shifty.==0.0)
  end

  @testset "Face data shift to dual cell center" begin
    cellx = Whirl.Grids.dualshiftx(dom.grid,facexunit)
    @test cellx[i,j] == 0.5 && cellx[i,j+1] == 0.5
    cellx[i,j] = cellx[i,j+1] = 0.0
    @test all(cellx.==0.0)

    celly = Whirl.Grids.dualshifty(dom.grid,faceyunit)
    @test celly[i,j] == 0.5 && celly[i+1,j] == 0.5
    celly[i,j] = celly[i+1,j] = 0.0
    @test all(celly.==0.0)
  end

  @testset "div curl" begin
    @test all(Whirl.Grids.diverg(dom.grid,Whirl.Grids.curl(dom.grid,cellunit)).==0.0)
  end

  @testset "curl grad" begin
    @test all(Whirl.Grids.curl(dom.grid,Whirl.Grids.grad(dom.grid,nodeunit)).==0.0)
  end

  Whirl.Grids.lgf_table!(dom.grid)

  @testset "Laplacian of the LGF" begin
    ψ = Whirl.Grids.L⁻¹(dom.grid)(cellunit);
    lapψ = Whirl.Grids.lap(dom.grid,ψ)
    @test lapψ[i,j]≈1.0
    @test isapprox(maximum(abs.(lapψ[Not(i),:])),0.0;atol=10.0*eps()) &&
            isapprox(maximum(abs.(lapψ[:,Not(j)])),0.0;atol=10.0*eps())
  end

end
