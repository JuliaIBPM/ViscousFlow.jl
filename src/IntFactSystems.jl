module IntFactSystems

export System, Constrained, Unconstrained
export A⁻¹, r₁

"Abstract type for an integrating factor system"
abstract type System{C} end

const Constrained = true
const Unconstrained = false

function A⁻¹ end
function r₁ end

end
