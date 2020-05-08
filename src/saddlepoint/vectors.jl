### Right-hand side and solution vectors

const SaddleVector = ArrayPartition

#SaddleVector(u::TU,f::TF) where {TU,TF} = ArrayPartition(u,f)

"""
    SaddleVector(u,f)

Construct a vector of a state part `u` and constraint part `f` of a
saddle-point vector, to be associated with a [`SaddleSystem`](@ref).
"""
function SaddleVector end

"""
    state(x::SaddleVector)

Provide the state part of the given saddle vector `x`
"""
state(u::SaddleVector) = u.x[1]

"""
    constraint(x::SaddleVector)

Provide the constraint part of the given saddle vector `x`
"""
constraint(u::SaddleVector) = u.x[2]
