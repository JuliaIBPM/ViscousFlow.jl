
struct PointForce
    x :: Tuple{Float64,Float64}
    f0 :: Float64
    t0 :: Float64
    σ :: Float64
    fbuffer :: ScalarData
    wbuffer :: Nodes
    regop :: Regularize
end

"""
    PointForce(x0::Tuple{Float64,Float64},f0,t0,σ,sys::NavierStokes)

Constructor function that immerses a point force in the dual nodes of system `sys`,
of strength `f0` to be applied at physical position `x0`, modulated by a Gaussian
centered at time `t0` with standard deviation `σ`.

The resulting function is a function of time and generates a field on dual nodes.
"""
function PointForce(x0::Tuple{Float64,Float64},f0::Float64,t0::Float64,σ::Float64,sys::NavierStokes)
    Xf = VectorData([x0[1]],[x0[2]])
    Ff = ScalarData(Xf)
    Ff[1] = f0
    wf = Nodes(Dual,size(sys))

    regop = Regularize(Xf,sys.Δx;I0=origin(sys),issymmetric=true)
    PointForce(x0,f0,t0,σ,Ff,wf,regop)
end

(f::PointForce)(t) = f.regop(f.wbuffer,rmul!(deepcopy(f.fbuffer),exp(-(t-f.t0)^2/f.σ^2)))
