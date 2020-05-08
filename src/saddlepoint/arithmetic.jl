### ARITHMETIC OPERATIONS

function mul!(output::Union{Tuple{AbstractVector{T},AbstractVector{T}},AbstractVectorOfArray{T}},sys::SaddleSystem{T,Ns,Nc},
  input::Union{Tuple{AbstractVector{T},AbstractVector{T}},AbstractVectorOfArray{T}}) where {T,Ns,Nc}

    @unpack A, B₂, B₁ᵀ, C = sys

    u,f = input
    r₁,r₂ = output
    length(u) == length(r₁) == Ns || error("Incompatible number of elements")
    length(f) == length(r₂) == Nc || error("Incompatible number of elements")

    r₁ .= A*u .+ B₁ᵀ*f
    r₂ .= B₂*u .+ C*f
    return output
end

function mul!(sol::Tuple{TU,TF},sys::SaddleSystem,rhs::Tuple{TU,TF}) where {T,Ns,Nc,TU,TF}
    u, f = sol
    r₁, r₂ = rhs
    return mul!((_unwrap_vec(u),_unwrap_vec(f)),sys,(_unwrap_vec(r₁),_unwrap_vec(r₂)))
end

function mul!(sol::ArrayPartition,sys::SaddleSystem,rhs::ArrayPartition) where {T,Ns,Nc,TU,TF}
    u, f = sol.x
    r₁, r₂ = rhs.x
    return mul!((_unwrap_vec(u),_unwrap_vec(f)),sys,(_unwrap_vec(r₁),_unwrap_vec(r₂)))
end

function (*)(sys::SaddleSystem,input::Tuple)
    u, f = input
    output = (similar(u),similar(f))
    mul!(output,sys,input)
    return output
end

# Routine for accepting vector inputs, parsing it into Ns and Nc parts
function mul!(sol::AbstractVector{T},sys::SaddleSystem{T,Ns,Nc},rhs::AbstractVector{T}) where {T,Ns,Nc}
    mul!(_split_vector(sol,Ns,Nc),sys,_split_vector(rhs,Ns,Nc))
    return sol
end

function (*)(sys::SaddleSystem{T},rhs::Union{AbstractVector{T},AbstractVectorOfArray{T},ArrayPartition{T}}) where {T}
    output = similar(rhs)
    mul!(output,sys,rhs)
    return output
end

## Multiplying tuples of saddle point systems

function mul!(sol,sys::NTuple{M,SaddleSystem},rhs) where {M}
   for (i,sysi) in enumerate(sys)
     mul!(sol[i],sysi,rhs[i])
   end
   sol
end

function (*)(sys::NTuple{M,SaddleSystem},rhs) where {M}
    sol = deepcopy.(rhs)
    mul!(sol,sys,rhs)
    return sol
end

#### Left division ####

function ldiv!(sol::Union{Tuple{AbstractVector{T},AbstractVector{T}},AbstractVectorOfArray{T}},sys::SaddleSystem{T,Ns,Nc},
              rhs::Union{Tuple{AbstractVector{T},AbstractVector{T}},AbstractVectorOfArray{T}}) where {T,Ns,Nc}
    @unpack A⁻¹, B₂, B₁ᵀ, B₂A⁻¹r₁, P, S⁻¹, _f_buf, A⁻¹B₁ᵀf = sys


    N = Ns+Nc
    u,f = sol
    r₁,r₂ = rhs
    length(u) == length(r₁) == Ns || error("Incompatible number of elements")
    length(f) == length(r₂) == Nc || error("Incompatible number of elements")

    u .= A⁻¹*r₁

    B₂A⁻¹r₁ .= B₂*u
    _f_buf .= r₂
    _f_buf .-= B₂A⁻¹r₁

    if Nc > 0
        f .= S⁻¹*_f_buf
        f .= P*f
    end
    A⁻¹B₁ᵀf .= A⁻¹*B₁ᵀ*f
    u .-= A⁻¹B₁ᵀf

    return sol
end

function ldiv!(sol::Tuple{TU,TF},sys::SaddleSystem,rhs::Tuple{TU,TF}) where {T,Ns,Nc,TU,TF}
    u, f = sol
    r₁, r₂ = rhs
    return ldiv!((_unwrap_vec(u),_unwrap_vec(f)),sys,(_unwrap_vec(r₁),_unwrap_vec(r₂)))
end

function ldiv!(sol::ArrayPartition,sys::SaddleSystem,rhs::ArrayPartition) where {T,Ns,Nc,TU,TF}
    u, f = sol.x
    r₁, r₂ = rhs.x
    return ldiv!((_unwrap_vec(u),_unwrap_vec(f)),sys,(_unwrap_vec(r₁),_unwrap_vec(r₂)))
end

function (\)(sys::SaddleSystem,rhs::Tuple) where {T,Ns,Nc}
    u, f = rhs
    sol = (similar(u),similar(f))
    ldiv!(sol,sys,rhs)
    return sol
end

# Routine for accepting vector inputs, parsing it into Ns and Nc parts
function ldiv!(sol::AbstractVector{T},sys::SaddleSystem{T,Ns,Nc},rhs::AbstractVector{T}) where {T,Ns,Nc}
    ldiv!(_split_vector(sol,Ns,Nc),sys,_split_vector(rhs,Ns,Nc))
    return sol
end

function (\)(sys::SaddleSystem{T},rhs::Union{AbstractVector{T},AbstractVectorOfArray{T},ArrayPartition{T}}) where {T}
    sol = similar(rhs)
    ldiv!(sol,sys,rhs)
    return sol
end

## Solving tuples of saddle point systems

function ldiv!(sol,sys::NTuple{M,SaddleSystem},rhs) where {M}
   for (i,sysi) in enumerate(sys)
     ldiv!(sol[i],sysi,rhs[i])
   end
   sol
end

function (\)(sys::NTuple{M,SaddleSystem},rhs) where {M}
    sol = deepcopy.(rhs)
    ldiv!(sol,sys,rhs)
    return sol
end

# vector -> tuple
_split_vector(x,Ns,Nc) = view(x,1:Ns), view(x,Ns+1:Ns+Nc)
