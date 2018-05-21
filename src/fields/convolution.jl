import Base: A_mul_B!, *

"""
    CircularConvolution{M, N}

A preplanned, circular convolution operator on an M × N matrix.

# Fields
- `Ĝ`: DFT coefficients of the convolution kernel
- `F`: preplanned rFFT operator
- `F⁻¹`: preplanned irFFT operator
- `paddedSpace`: scratch space to zero-pad the input matrix
- `Â`: scratch space to store the DFT coefficients of the zero-padded input matrix

# Constructors:

- `CircularConvolution(G::Matrix{Float64})`

# Example:
```jldoctest
julia> G = repmat(1.0:3,1,4)
3×4 Array{Int64,2}:
 1  1  1  1
 2  2  2  2
 3  3  3  3

julia> C = CircularConvolution(G)
Zero-padded convolution for a 3 × 4 matrix

julia> C*reshape(1:12, 3, 4)
3×4 Array{Int64,2}:
 164  164  164  164
 130  130  130  130
 148  148  148  148
```
"""
struct CircularConvolution{M, N, K, KI}
    Ĝ::Matrix{Complex128}
    F::K
    F⁻¹::KI

    paddedSpace::Matrix{Float64}
    Â::Matrix{Complex128}
end

function Base.show(io::IO, c::CircularConvolution{M, N}) where {M, N}
    print(io, "Circular convolution on a $M × $N matrix")
end

function CircularConvolution(G::AbstractMatrix{Float64}, fftw_flags = FFTW.ESTIMATE)
    M, N = size(G)
    paddedSpace = Matrix{Float64}(2M-1, 2N-1)
    F = FFTW.plan_rfft(paddedSpace, flags = fftw_flags)

    mirror!(paddedSpace, G)
    Ĝ = F * paddedSpace

    Â = similar(Ĝ)
    F⁻¹ = FFTW.plan_irfft(Â, 2M - 1, flags = fftw_flags)

    CircularConvolution{M, N, typeof(F), typeof(F⁻¹)}(Ĝ, F, F⁻¹, paddedSpace, Â)
end

function A_mul_B!(out, C::CircularConvolution{M, N}, B) where {M, N}
    @assert size(out) == size(B) == (M, N)

    inds = CartesianRange((M,N))
    fill!(C.paddedSpace, 0)
    copy!(C.paddedSpace, inds, B, inds)
    A_mul_B!(C.Â, C.F, C.paddedSpace)

    C.Â .*= C.Ĝ

    A_mul_B!(C.paddedSpace, C.F⁻¹, C.Â)

    copy!(out, inds, C.paddedSpace, CartesianRange((M:2M-1,N:2N-1)))
end

C::CircularConvolution * B = A_mul_B!(similar(B), C, B)

function mirror!(A, a::AbstractArray{T,2}) where {T}
    Nr, Nc = size(a)
    A[1:Nr-1, 1:Nc-1] .= a[Nr:-1:2, Nc:-1:2]
    A[1:Nr-1, Nc:end] .= a[Nr:-1:2, 1:Nc]
    A[Nr:end, 1:Nc-1] .= a[1:Nr, Nc:-1:2]
    A[Nr:end, Nc:end] .= a
    A
end

function mirror(a::AbstractArray{T,2}) where {T}
    Nr, Nc = size(a)
    mirror!(zeros(T, 2Nr-1, 2Nc-1), a)
end
