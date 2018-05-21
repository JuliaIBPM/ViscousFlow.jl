module IntFactSystems

export System, Constrained, Unconstrained
export A⁻¹, r₁, intfact

"Abstract type for an integrating factor system"
abstract type System{C} end

const Constrained = true
const Unconstrained = false

function A⁻¹ end
function r₁ end

intfact(x, y,a) = exp(-4a)besseli(x,2a)besseli(y,2a)

end
