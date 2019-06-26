
# For systems of rigid bodies, this constructs the right-hand side of the equations
# with the velocity prescribed by the motions in ml
function TimeMarching.r₁(u::NTuple{NB,Vector{Float64}},t::Real,ml::Vector{RigidBodyMotion}) where {NB}
    du = deepcopy(u)
    for ib = 1:NB
        _,ċ,_,_,α̇,_ = ml[ib](t)
        du[ib] .= [real(ċ),imag(ċ),α̇]
    end
    return du
end
