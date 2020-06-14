using LinearAlgebra

@testset "PointForce " begin

  @testset "Gaussian function" begin
    @test Gaussian!(1.0, 0.1) ≈ exp(-0.5*1.0^2*0.1^(-2))
    @test Gaussian(2.0)(0.5) ≈ Gaussian!(0.5, 2.0)
  end

  @testset "Spatial Gaussian distribution for NodeData" begin

    Re = 500; # Reynolds number
    xlim = (-1.0,5.0)
    ylim = (-2.0,2.0)
    ReΔ = 16
    Δx = ReΔ/Re
    Δt = min(0.25*Δx,0.5*Δx^2*Re);
    sys = NavierStokes(Re, Δx, xlim, ylim, Δt)

    w = Nodes(Dual,size(sys));
    wbuffer =  Nodes(Dual,size(sys));
    x0 = (1.0, 2.0)
    σx = 2.0

    SpatialGauss(wbuffer, w, x0, σx, sys)
    xx, yy = coordinates(w,dx=sys.grid.Δx,I0=origin(sys))
    @test wbuffer[10, 20] ≈ exp(-0.5*σx^(-2)*((xx[10]-x0[1])^2+(yy[20]-x0[2])^2))
  end

  @testset "Spatial Gaussian distribution for EdgeData" begin

    Re = 500; # Reynolds number
    xlim = (-1.0,5.0)
    ylim = (-2.0,2.0)
    ReΔ = 16
    Δx = ReΔ/Re
    Δt = min(0.25*Δx,0.5*Δx^2*Re);
    sys = NavierStokes(Re, Δx, xlim, ylim, Δt)

    q = Edges(Primal,size(sys));
    qbuffer =  Edges(Primal,size(sys));
    x0 = (1.0, 2.0)
    σx = 2.0

    SpatialGauss(qbuffer, q, x0, σx, sys)
    xx, xy, yx, yy = coordinates(q,dx=sys.grid.Δx,I0=origin(sys))
    @test qbuffer.u[10, 20] ≈ exp(-0.5*σx^(-2)*((xx[10]-x0[1])^2+(xy[20]-x0[2])^2))

    @test qbuffer.v[60, 41] ≈ exp(-0.5*σx^(-2)*((yx[60]-x0[1])^2+(yy[41]-x0[2])^2))

  end

@testset "PointForce(t) without spatial variation" begin
  Re = 500; # Reynolds number
  xlim = (-1.0,5.0)
  ylim = (-2.0,2.0)
  ReΔ = 16
  Δx = ReΔ/Re
  Δt = min(0.25*Δx,0.5*Δx^2*Re);
  sys = NavierStokes(Re, Δx, xlim, ylim, Δt)

  q = Edges(Primal,size(sys));
  qbuffer =  Edges(Primal,size(sys));
  xx, xy, yx, yy = coordinates(q,dx=sys.grid.Δx,I0=origin(sys))
  x0 = (xx[10], xy[30])
  f0 = (1.0, 2.0)
  σx = 2.0
  t0 = 2.0
  σt = 0.5
  pointf = PointForce(q, x0, f0, t0, σt, sys)

  # Test at the central location and time of the Gaussian distribution
  @test pointf(2.0).u[10, 30] ≈ (pointf.regop(pointf.ubuffer,deepcopy(pointf.fbuffer))).u[10, 30]

  # Test at the central time of the Gaussian distribution
  @test pointf(2.0).v.data  ≈  pointf.regop(pointf.ubuffer,deepcopy(pointf.fbuffer)).v.data

  @test pointf(2.5).u.data  ≈  pointf.regop(pointf.ubuffer,rmul!(deepcopy(pointf.fbuffer),
                              exp(-(2.5-pointf.t0)^2/pointf.σt^2))).u.data

  @test pointf(2.1).v.data  ≈  pointf.regop(pointf.ubuffer,rmul!(deepcopy(pointf.fbuffer),
                              exp(-(2.1-pointf.t0)^2/pointf.σt^2))).v.data

end

@testset "PointForce(t) with spatial variation" begin
        Re = 500; # Reynolds number
        xlim = (-1.0,5.0)
        ylim = (-2.0,2.0)
        ReΔ = 16
        Δx = ReΔ/Re
        Δt = min(0.25*Δx,0.5*Δx^2*Re);
        sys = NavierStokes(Re, Δx, xlim, ylim, Δt)

        q = Edges(Primal,size(sys));
        qbuffer =  Edges(Primal,size(sys));
        xx, xy, yx, yy = coordinates(q,dx=sys.grid.Δx,I0=origin(sys))
        x0 = (xx[10], xy[30])
        f0 = (1.0, 2.0)
        σx = 2.0
        t0 = 2.0
        σt = 0.5
        pointf = PointForce(q, x0, f0, t0, σt, sys, σx = σx)

        # Test at the central location and time of the Gaussian distribution
        @test pointf(2.0).u[10, 30] ≈  (pointf.regop(pointf.ubuffer,deepcopy(pointf.fbuffer))).u[10, 30]

        # Test at the central time of the Gaussian distribution
        @test pointf(2.0).v.data ≈  (pointf.xbuffer∘pointf.regop(pointf.ubuffer,deepcopy(pointf.fbuffer))).v.data

        @test pointf(2.5).u.data ≈  (pointf.xbuffer∘pointf.regop(pointf.ubuffer,rmul!(deepcopy(pointf.fbuffer),
                                    exp(-(2.5-pointf.t0)^2/pointf.σt^2)))).u.data

        @test pointf(2.1).v.data ≈  (pointf.xbuffer∘pointf.regop(pointf.ubuffer,rmul!(deepcopy(pointf.fbuffer),
                                    exp(-(2.1-pointf.t0)^2/pointf.σt^2)))).v.data

  end
end
