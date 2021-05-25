#=
Create caches for terms in the equations here.
=#

abstract type AbstractPDECache end

struct ConvectiveDerivativeCache{VT,DVT} <: AbstractPDECache
  V :: VT
  Vtf :: DVT
  DVf :: DVT
  VDVf :: DVT
end

struct RHSCache{VT} <: AbstractPDECache
  dv :: VT
end

struct PressureCache{VT,PT} <: AbstractPDECache
  dvp :: VT
  dvf :: VT
  lp :: PT
end

struct VelocityCache{WT,VT,FT,VBT,SBT} <: AbstractPDECache
  Wn :: WT
  Vf :: VT
  Sc :: FT
  Δus :: VBT
  Sb :: SBT
end

struct StreamfunctionCache{ST} <: AbstractPDECache
  Sn :: ST
end

struct ConstraintOperatorCache{VT,VBT} <: AbstractPDECache
  Vv :: VT
  Δus :: VBT
end

struct DoubleLayerCache{VT,VBT} <: AbstractPDECache
  Vv :: VT
  Δus :: VBT
end

struct SurfaceVelocityCache{VBT} <: AbstractPDECache
  Vb :: VBT
end

struct SurfaceTractionCache{VBT,SBT} <: AbstractPDECache
  Vb :: VBT
  Δus :: VBT
  Sb :: SBT
end
