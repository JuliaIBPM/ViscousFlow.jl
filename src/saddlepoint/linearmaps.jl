### LINEAR MAP CONSTRUCTION

# for a given function of function-like object A, which acts upon data of type `input`
# and returns data of type `output`
# return a LinearMap that acts upon a vector form of u
linear_map(A,input,output;eltype=Float64) = _linear_map(A,input,output,eltype,Val(length(input)),Val(length(output)))

linear_map(A,input;eltype=Float64) = _linear_map(A,input,eltype,Val(length(input)))

linear_map(A::AbstractMatrix{T},input::AbstractVector{T};eltype=Float64) where {T} = LinearMap{eltype}(A)

linear_map(A::SaddleSystem{T,Ns,Nc},::Any;eltype=Float64) where {T,Ns,Nc} = LinearMap{T}(x->A*x,Ns+Nc)

linear_map(A::UniformScaling,input;eltype=Float64) = LinearMap{eltype}(x->A*x,length(input))

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
      LinearMap{eltype}(x -> _unwrap_vec(zero(output)),length(output),0)

# output is 0 length, input is not
_linear_map(A,input,output,eltype,::Val{N},::Val{0}) where {N} =
      LinearMap{eltype}(x -> (),0,length(input))

# non-zero lengths of input and output
_linear_map(A,input,output,eltype,::Val{N},::Val{M}) where {N,M} =
      LinearMap{eltype}(_create_fcn(A,input),length(output),length(input))


# Create a function for operator A that can act upon an input of type AbstractVector
# and return an output of type AbstractVector. It should wrap the input vector
# in the input data type associated with A and it should then unwrap its
# output back into vector form
function _create_fcn(A,input)
    # if A has an associated * operation, then use this
    if hasmethod(*,Tuple{typeof(A),typeof(input)})
        fcn = _create_vec_multiplication(A,input)
    # or if A is a function or function-like object, then use this
    elseif hasmethod(A,Tuple{typeof(input)})
        fcn = _create_vec_function(A,input)
    # or just quit
    else
        error("No function exists for this operator to act upon this type of data")
    end
    return fcn
end

_create_vec_multiplication(A,u::TU) where {TU} = (x -> _unwrap_vec(A*_wrap_vec(x,u)))
_create_vec_function(A,u::TU) where {TU} = (x -> _unwrap_vec(A(_wrap_vec(x,u))))
_create_vec_backslash(A,u::TU) where {TU} = (x -> _unwrap_vec(A\_wrap_vec(x,u)))

#### WRAPPERS ####
# wrap the vector x in type u, unless u is already a subtype of AbstractVector,
# in which case it just keeps it as is.
_wrap_vec(x::AbstractVector{T},u::TU) where {T,TU} = TU(reshape(x,size(u)...))

# if the vector x is simply a reshaped form of type u, then just get the
# parent of x
_wrap_vec(x::Base.ReshapedArray,u::TU) where {TU} = parent(x)


#### UNWRAPPERS ####
# Usually vec suffices
_unwrap_vec(x) = vec(x)
_unwrap_vec(x::Tuple) = x
