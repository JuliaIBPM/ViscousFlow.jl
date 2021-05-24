#=
Create subtypes for all terms in the equations here.
The goal of this should be to change NavierStokes into
a more general type for any PDEs, similar to the structures
used for OrdinaryDiffEq.

Rule of thumb: keep the typing generic!

Other things to keep in mind:
- For MovingPoints problems, the immersion operators will get changed. These
   need to be easily accessible from the overall system. We can endow
   the term type with a parameter for types that must be updated if the points
   update. Each such type will need to be accompanied by an operator for
   setting it up (and possibly a separate one for updating it).
- Some operators, such as Rf and Ef, may get used by more than one term.
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

struct PressureCache{DVT,PT} <: AbstractPDECache
  dvp :: DVT
  dvf :: DVT
  lp :: PT
end
