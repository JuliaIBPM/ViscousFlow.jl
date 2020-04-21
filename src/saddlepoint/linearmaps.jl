### LINEAR MAP CONSTRUCTION

# for a given function of function-like object A, which acts upon data of type u
# and returns data of type f
# return a LinearMap that acts upon a vector form of u
linear_map(A,u,f;eltype=Float64) = _linear_map(A,u,f,eltype)
linear_map(A,u;eltype=Float64) = _linear_map(A,u,eltype)
linear_map(A::AbstractMatrix{T},u::AbstractVector{T};eltype=Float64) where {T} = LinearMap{eltype}(A)
linear_map(A::SaddleSystem{T,Ns,Nc},::Any;eltype=Float64) where {T,Ns,Nc} = LinearMap{T}(x->A*x,Ns+Nc)

function linear_inverse_map(A,input;eltype=Float64)
    hasmethod(\,Tuple{typeof(A),typeof(input)}) || error("No such backslash operator exists")
    return LinearMap{eltype}(_create_vec_backslash(A,input),length(input))
end
linear_inverse_map(A::SaddleSystem{T,Ns,Nc},::Any;eltype=Float64) where {T,Ns,Nc} = LinearMap{T}(x->A\x,Ns+Nc)

_linear_map(A,input,output,eltype) = LinearMap{eltype}(_create_fcn(A,input),length(output),length(input))
_linear_map(A,input,eltype) = LinearMap{eltype}(_create_fcn(A,input),length(input))

function _create_fcn(A,input)
    if hasmethod(*,Tuple{typeof(A),typeof(input)})
        fcn = _create_vec_multiplication(A,input)
    elseif hasmethod(A,Tuple{typeof(input)})
        fcn = _create_vec_function(A,input)
    end
    return fcn
end

_create_vec_multiplication(A,u::TU) where {TU} = (x -> vec(A*_wrap_vec(x,u)))
_create_vec_function(A,u::TU) where {TU} = (x -> vec(A(_wrap_vec(x,u))))
_create_vec_backslash(A,u::TU) where {TU} = (x -> vec(A\_wrap_vec(x,u)))

_wrap_vec(x::Vector{T},u::TU) where {T,TU} = TU(reshape(x,size(u)...))
_wrap_vec(x::Base.ReshapedArray,u::TU) where {TU} = parent(x)
_wrap_vec(x::AbstractVector{T},u::TU) where {T,TU} = x
_wrap_vec(x,u::TU) where {TU <: Tuple} = x

_unwrap_vec(x) = vec(x)
_unwrap_vec(x::Tuple) = x
