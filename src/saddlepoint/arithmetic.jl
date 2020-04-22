### ARITHMETIC OPERATIONS

function mul!(output::Tuple{AbstractVector{T},AbstractVector{T}},sys::SaddleSystem{T,Ns,Nc},input::Tuple{AbstractVector{T},AbstractVector{T}}) where {T,Ns,Nc}
    u,f = input
    r₁,r₂ = output
    length(u) == length(r₁) == Ns || error("Incompatible number of elements")
    length(f) == length(r₂) == Nc || error("Incompatible number of elements")

    r₁ .= sys.A*u + sys.B₁ᵀ*f
    r₂ .= sys.B₂*u + sys.C*f
    return output
end

function mul!(sol::Tuple{TU,TF},sys::SaddleSystem,rhs::Tuple{TU,TF}) where {T,Ns,Nc,TU,TF}
    u, f = sol
    r₁, r₂ = rhs
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

function (*)(sys::SaddleSystem,input::AbstractVector)
    output = similar(input)
    mul!(output,sys,input)
    return output
end

# Left division

function ldiv!(sol::Tuple{AbstractVector{T},AbstractVector{T}},sys::SaddleSystem{T,Ns,Nc},rhs::Tuple{AbstractVector{T},AbstractVector{T}}) where {T,Ns,Nc}

    N = Ns+Nc
    u,f = sol #_split_vector(sol,Ns,Nc)
    r₁,r₂ = rhs #_split_vector(rhs,Ns,Nc)
    length(u) == length(r₁) == Ns || error("Incompatible number of elements")
    length(f) == length(r₂) == Nc || error("Incompatible number of elements")

    u .= sys.A⁻¹*r₁

    sys.B₂A⁻¹r₁ .= sys.B₂*u
    sys._f_buf .= r₂
    sys._f_buf .-= sys.B₂A⁻¹r₁

    if Nc > 0
        f .= sys.S⁻¹*sys._f_buf
        #f .= sys.P(f)
    end
    sys.A⁻¹B₁ᵀf .= sys.A⁻¹*sys.B₁ᵀ*f
    u .-= sys.A⁻¹B₁ᵀf

    return sol
end

function ldiv!(sol::Tuple{TU,TF},sys::SaddleSystem,rhs::Tuple{TU,TF}) where {T,Ns,Nc,TU,TF}
    u, f = sol
    r₁, r₂ = rhs
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

function (\)(sys::SaddleSystem,rhs::AbstractVector)
    sol = similar(rhs)
    ldiv!(sol,sys,rhs)
    return sol
end

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
