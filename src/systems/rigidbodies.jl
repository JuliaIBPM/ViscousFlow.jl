
# For systems of rigid bodies, this constructs the right-hand side of the equations
# with the velocity prescribed by the motions in ml
# The state vector is a concatenation of the 3 x 1 vectors for each body
function TimeMarching.r₁(u::Vector{Float64},t::Real,ml::Vector{RigidBodyMotion})
    du = deepcopy(u)
    cnt = 0
    for ib = 1:length(ml)
        _,ċ,_,_,α̇,_ = ml[ib](t)
        du[cnt+1:cnt+3] = [real(ċ),imag(ċ),α̇]
        cnt += 3
    end
    return du
end
