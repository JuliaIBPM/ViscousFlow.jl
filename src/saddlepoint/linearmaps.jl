### LINEAR MAP CONSTRUCTION

# for a given function of function-like object A, which acts upon data of type `input`
# and returns data of type `output`
# return a LinearMap that acts upon a vector form of u
linear_map(A,input,output;eltype=Float64) = _linear_map(A,input,output,eltype,Val(length(input)),Val(length(output)))

linear_map(A,input;eltype=Float64) = _linear_map(A,input,eltype,Val(length(input)))

linear_map(A::AbstractMatrix{T},input::AbstractVector{T};eltype=Float64) where {T} = LinearMap{eltype}(A)

linear_map(A::SaddleSystem{T,Ns,Nc},::Any;eltype=Float64) where {T,Ns,Nc} = LinearMap{T}(x->A*x,Ns+Nc)

function linear_inverse_map(A,input;eltype=Float64)
    hasmethod(\,Tuple{typeof(A),typeof(input)}) || error("No such backslash operator exists")
    return LinearMap{eltype}(_create_vec_backslash(A,input),length(input))
end
linear_inverse_map(A::SaddleSystem{T,Ns,Nc},::Any;eltype=Float64) where {T,Ns,Nc} = LinearMap{T}(x->A\x,Ns+Nc)

# Square operators. input of zero length
_linear_map(A,input,eltype,::Val{0}) =
      LinearMap{eltype}(x -> (),0,0)

# Square operators. input of non-zero length
_linear_map(A,input,eltype,::Val{M}) where {M} =
      LinearMap{eltype}(_create_fcn(A,input),length(input))

# input and output have zero lengths
_linear_map(A,input,output,eltype,::Val{0},::Val{0}) =
      LinearMap{eltype}(x -> (),0,0)

# input is 0 length, output is not
_linear_map(A,input,output,eltype,::Val{0},::Val{M}) where {M} =
      LinearMap{eltype}(x -> _unwrap_vec(typeof(output)()),length(output),0)

# output is 0 length, input is not
_linear_map(A,input,output,eltype,::Val{N},::Val{0}) where {N} =
      LinearMap{eltype}(x -> (),0,length(input))

# non-zero lengths of input and output
_linear_map(A,input,output,eltype,::Val{N},::Val{M}) where {N,M} =
      LinearMap{eltype}(_create_fcn(A,input),length(output),length(input))




function _create_fcn(A,input)
    if hasmethod(*,Tuple{typeof(A),typeof(input)})
        fcn = _create_vec_multiplication(A,input)
    elseif hasmethod(A,Tuple{typeof(input)})
        fcn = _create_vec_function(A,input)
    end
    return fcn
end

_create_vec_multiplication(A,u::TU) where {TU} = (x -> _unwrap_vec(A*_wrap_vec(x,u)))
_create_vec_function(A,u::TU) where {TU} = (x -> _unwrap_vec(A(_wrap_vec(x,u))))
_create_vec_backslash(A,u::TU) where {TU} = (x -> _unwrap_vec(A\_wrap_vec(x,u)))

# wrap the vector x in type u, unless u is already a subtype of AbstractVector
_wrap_vec(x::AbstractVector{T},u::TU) where {T,TU} = TU(reshape(x,size(u)...))
_wrap_vec(x::AbstractVector{T},u::TU) where {T,TU <: AbstractVector} = x

# if the vector x is simply a reshaped form of type u, then just get the
# parent of x
_wrap_vec(x::Base.ReshapedArray,u::TU) where {TU} = parent(x)

# not sure if this one is needed
#_wrap_vec(x,u::TU) where {TU <: Tuple} = x

_unwrap_vec(x) = vec(x)
_unwrap_vec(x::Tuple) = x
