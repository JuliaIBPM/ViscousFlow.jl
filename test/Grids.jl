@testset "Grid Routines" begin

  xmin = [-1.0,-1.0]
  xmax = [1.0,1.0]
  dom = Whirl.Systems.DualDomain(xmin,xmax)

  Δx = 0.2
  dom = Whirl.Systems.add_grid(dom,Δx)

  i = 5; j = 5

  cellzero = zeros(dom.grid.cell)
  nodezero = zeros(dom.grid.node)
  facexzero = zeros(dom.grid.facex)
  faceyzero = zeros(dom.grid.facey)
  dualfacexzero = zeros(dom.grid.dualfacex)
  dualfaceyzero = zeros(dom.grid.dualfacex)

  @testset "Dual cell center data Laplacian" begin
    cellunit = copy(cellzero)
    cellunit[i,j] = 1.0
    lapcell = Whirl.Grids.lap(dom.grid,cellunit)
    @test lapcell[i,j] == -4.0
    @test lapcell[i+1,j] == lapcell[i-1,j] == lapcell[i,j-1] == lapcell[i,j+1] == 1.0
  end

  @testset "Dual cell center data curl" begin
    cellunit = copy(cellzero)
    cellunit[i,j] = 1.0
    qx,qy = Whirl.Grids.curl(dom.grid,cellunit)
    @test qx[i,j-1] == 1.0 && qx[i,j] == -1.0
    @test qy[i-1,j] == -1.0 && qy[i,j] == 1.0
  end

  @testset "Dual cell node gradient" begin
    nodeunit = copy(nodezero)
    nodeunit[i,j] = 1.0
    qx,qy = Whirl.Grids.grad(dom.grid,nodeunit)
    @test qx[i,j] == 1.0 && qx[i+1,j] == -1.0
    @test qy[i,j] == 1.0 && qy[i,j+1] == -1.0
  end

  @testset "Face data curl" begin
    facexunit = copy(facexzero)
    facexunit[i,j] = 1.0
    cellcurl = Whirl.Grids.curl(dom.grid,facexunit,faceyzero)
    @test cellcurl[i,j] == -1.0 && cellcurl[i,j+1] == 1.0
    faceyunit = copy(faceyzero)
    faceyunit[i,j] = 1.0
    cellcurl = Whirl.Grids.curl(dom.grid,facexzero,faceyunit)
    @test cellcurl[i,j] == 1.0 && cellcurl[i+1,j] == -1.0
  end

  @testset "Face data divergence" begin
    facexunit = copy(facexzero)
    facexunit[i,j] = 1.0
    nodediv = Whirl.Grids.diverg(dom.grid,facexunit,faceyzero)
    @test nodediv[i,j] == -1.0 && nodediv[i-1,j] == 1.0
    faceyunit = copy(faceyzero)
    faceyunit[i,j] = 1.0
    nodediv = Whirl.Grids.diverg(dom.grid,facexzero,faceyunit)
    @test nodediv[i,j] == -1.0 && nodediv[i,j-1] == 1.0
  end

  @testset "Face data Laplacian" begin
    facexunit = copy(facexzero)
    facexunit[i,j] = 1.0
    lapx,lapy = Whirl.Grids.lap(dom.grid,facexunit,faceyzero)
    @test lapx[i,j] == -4.0
    @test lapx[i+1,j] == lapx[i-1,j] == lapx[i,j-1] == lapx[i,j+1] == 1.0
    @test maximum(abs.(lapy)) == 0.0

    faceyunit = copy(faceyzero)
    faceyunit[i,j] = 1.0
    lapx,lapy = Whirl.Grids.lap(dom.grid,facexzero,faceyunit)
    @test lapy[i,j] == -4.0
    @test lapy[i+1,j] == lapy[i-1,j] == lapy[i,j-1] == lapy[i,j+1] == 1.0
    @test maximum(abs.(lapx)) == 0.0
  end

  @testset "Dual face data divergence" begin
    dualfacexunit = copy(dualfacexzero)
    dualfacexunit[i,j] = 1.0
    celldiv = Whirl.Grids.dualdiverg(dom.grid,dualfacexunit,dualfaceyzero)
    @test celldiv[i,j] == 1.0 && celldiv[i+1,j] == -1.0

    dualfaceyunit = copy(dualfaceyzero)
    dualfaceyunit[i,j] = 1.0
    celldiv = Whirl.Grids.dualdiverg(dom.grid,dualfacexzero,dualfaceyunit)
    @test celldiv[i,j] == 1.0 && celldiv[i,j+1] == -1.0
  end

  @testset "Face data shift to dual face" begin
    facexunit = copy(facexzero)
    facexunit[i,j] = 1.0
    shiftx,shifty = Whirl.Grids.shift(dom.grid,facexunit,faceyzero)
    @test shiftx[i,j] == shiftx[i-1,j] == shiftx[i,j+1] == shiftx[i-1,j+1] == 0.25
  end

  @testset "Dual cell center data shift to dual face" begin
    cellunit = copy(cellzero)
    cellunit[i,j] = 1.0
    shiftx,shifty = Whirl.Grids.shift(dom.grid,cellunit)
    @test shiftx[i,j] == shiftx[i-1,j] == 0.5
    @test shifty[i,j] == shifty[i,j-1] == 0.5
  end

  @testset "Face data shift to dual cell center" begin
    facexunit = copy(facexzero)
    facexunit[i,j] = 1.0
    cellx = Whirl.Grids.dualshiftx(dom.grid,facexunit)
    @test cellx[i,j] == 0.5 && cellx[i,j+1] == 0.5

    faceyunit = copy(faceyzero)
    faceyunit[i,j] = 1.0
    celly = Whirl.Grids.dualshifty(dom.grid,faceyunit)
    @test celly[i,j] == 0.5 && celly[i+1,j] == 0.5
  end



end
