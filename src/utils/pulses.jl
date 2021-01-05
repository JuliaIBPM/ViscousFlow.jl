### Utilities ###

struct PulseParams{GF}
    field :: GF
    t0 :: Float64
    σt :: Float64
end

_process_pulses(::Nothing,s,grid) = nothing

function _process_pulses(params::PulseParams,s,grid)
  gf = GeneratedField(s,params.field,grid)
  pulse = PulseField(gf,params.t0,params.σt)
  return [pulse]
end

function _process_pulses(params::Vector{<:PulseParams},s,grid)
  pulses = PulseField[]
  for p in params
    gf = GeneratedField(s,p.field,grid)
    push!(pulses,PulseField(gf,p.t0,p.σt))
  end
  return pulses
end

# The stuff below should be added to CartesianGrids eventually
import CartesianGrids: SpatialGaussian

struct DGaussian
  σ :: Float64
  x0 :: Float64
  A :: Float64
  fact :: Float64
end
DGaussian(σ,x0,A) = DGaussian(σ,x0,A, A/sqrt(π)/σ^2)

radius(g::DGaussian) = g.σ
center(g::DGaussian) = g.x0
strength(g::DGaussian) = g.A

@inline dgaussian(r;tol=6.0) = abs(r) < tol ? -2*r*exp(-r^2) : 0.0

(g::DGaussian)(x) = g.fact*dgaussian((x-center(g))/radius(g))

struct SpatialDGaussian{GX,GY} <: CartesianGrids.AbstractSpatialField
  gx :: GX
  gy :: GY
  A :: Float64
  SpatialDGaussian(gx,gy,A) = new{typeof(gx),typeof(gy)}(gx,gy,A)
end

"""
    SpatialGaussian(σx,σy,x0,y0,A[,derivdir=0])

Set up a spatial field in the form of a Gaussian centered at `x0,y0` with
radii `σx` and `σy` in the respective directions and amplitude `A`. If the
optional parameter `deriv` is set to 1 or 2, then it returns the first
derivative of a Gaussian in that direction (`x` or `y`, respectively).
"""
SpatialGaussian(σx,σy,x0,y0,A;deriv::Int=0) = _spatialdgaussian(σx,σy,x0,y0,A,Val(deriv))

_spatialdgaussian(σx,σy,x0,y0,A,::Val{0}) = SpatialGaussian(Gaussian(σx,x0,A),Gaussian(σy,y0,1),A)
_spatialdgaussian(σx,σy,x0,y0,A,::Val{1}) = SpatialDGaussian(DGaussian(σx,x0,A),Gaussian(σy,y0,1),A)
_spatialdgaussian(σx,σy,x0,y0,A,::Val{2}) = SpatialDGaussian(Gaussian(σx,x0,A),DGaussian(σy,y0,1),A)


SpatialDGaussian(σ,x0,y0,A,dir::Int) = SpatialDGaussian(σ,σ,x0,y0,A,dir)


(g::SpatialDGaussian)(x,y) = g.gx(x)*g.gy(y)
