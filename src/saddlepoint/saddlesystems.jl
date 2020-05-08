### SaddleSystem ###

abstract type SchurSolverType end

struct SaddleSystem{T,Ns,Nc,TF,TU,TS<:SchurSolverType}
    A :: LinearMap{T}
    B₂ :: LinearMap{T}
    B₁ᵀ :: LinearMap{T}
    C :: LinearMap{T}
    A⁻¹ :: LinearMap{T}
    A⁻¹B₁ᵀf :: Vector{T}
    B₂A⁻¹r₁ :: Vector{T}
    _f_buf :: Vector{T}
    P :: LinearMap{T}
    S :: LinearMap{T}
    S⁻¹ :: LinearMap{T}
end


# Constructors

"""
    SaddleSystem

Construct a saddle-point system operator from the constituent operator blocks. The resulting object can be used
with `*` and `\\` to multiply and solve. The saddle-point problem has the form

``
\\begin{bmatrix}A & B_1^T \\\\ B_2 & C \\end{bmatrix} \\begin{pmatrix} u \\\\ f \\end{pmatrix} = \\begin{pmatrix} r_1 \\\\ r_2 \\end{pmatrix}
``

### Constructors
`SaddleSystem(A::AbstractMatrix,B₂::AbstractMatrix,B₁ᵀ::AbstractMatrix,C::AbstractMatrix[,eltype=Float64])`.
Blocks are given as matrices. Must have consistent sizes to stack appropriately. If
this is called with `SaddleSystem(A,B₂,B₁ᵀ)`, it sets `C` to zero automatically.

`SaddleSystem(A,B₂,B₁ᵀ,C,u,f[,eltype=Float64])`.
Operators `A`, `B₂`, `B₁ᵀ`, `C` are given in various forms, including matrices, functions, and function-like objects.
`u` and `f` are examples of the data types in the corresponding solution and right-hand side vectors.
Guidelines:

* The entries `A` and `B₂` must be able to act upon `u` (either by multiplication or as a function) and `B₁ᵀ` and `C` must be able to act on `f` (also, either by multiplication or as a function).
* `A` and `B₁ᵀ` should return data of type `u`, and `B₂` and `C` should return data of type `f`.
* `A` must be invertible and be outfitted with operators `\` and `ldiv!`.
* Both `u` and `f` must be subtypes of `AbstractArray`: they must be equipped with `size`
  and `vec` functions and with a constructor of the form `T(data)` where `T` is the data type of
  `u` or `f` and `data` is the wrapped data array.

If called as `SaddleSystem(A,B₂,B₁ᵀ,u,f)`, the `C` block is omitted and assumed to be zero.

If called with `SaddleSystem(A,u)`, this is equivalent to calling `SaddleSystem(A,nothing,nothing,u,[])`, then this reverts
to the unconstrained system described by operator `A`.

The list of vectors `u` and `f` in any of these constructors can be bundled together
as a [`SaddleVector`](@ref), e.g. `SaddleSystem(A,B₂,B₁ᵀ,SaddleVector(u,f))`.
"""
function SaddleSystem(A::LinearMap{T},B₂::LinearMap{T},B₁ᵀ::LinearMap{T},C::LinearMap{T},
                      A⁻¹::LinearMap{T},P::LinearMap{T},TU,TF;solver::Type{TS}=Direct,kwargs...) where {T,TS<:SchurSolverType}

    ns, nc = _check_sizes(A,B₂,B₁ᵀ,C,P)

    S = C - B₂*A⁻¹*B₁ᵀ

    S⁻¹ = _schur_inverse_function(S,T,nc,solver,kwargs...)

    return SaddleSystem{T,ns,nc,TU,TF,solver}(A,B₂,B₁ᵀ,C,A⁻¹,zeros(T,ns),zeros(T,nc),zeros(T,nc),P,S,S⁻¹)
end

##### Schur complement solver functions #####

abstract type Direct <: SchurSolverType end

function _schur_inverse_function(S,T,M,::Type{Direct},kwargs...)
  Sfact = factorize(Matrix(S))
  return LinearMap{T}(x -> Sfact\x,M)
end

macro createsolver(stype)
  sroutine = Symbol(lowercase(string(stype)),"!")

  return esc(quote
          export $stype
          abstract type $stype <: SchurSolverType end
          function _schur_inverse_function(S,T,M,::Type{$stype},kwargs...)
            return LinearMap{T}(x -> (y = deepcopy(x); $sroutine(y,S,x;kwargs...); return y),M)
          end
        end)
end

@createsolver CG
@createsolver BiCGStabl
@createsolver GMRES
@createsolver MINRES
@createsolver IDRS

###########


### OTHER CONSTRUCTORS

### Matrix operators
function SaddleSystem(A::AbstractMatrix{T},B₂::AbstractMatrix{T},B₁ᵀ::AbstractMatrix{T},
                      C::AbstractMatrix{T};solver::Type{TS}=Direct,filter=I,kwargs...) where {T,TS<:SchurSolverType}

    Afact = factorize(A)
    Ainv = LinearMap{T}(x -> Afact\x,size(A,1))

    return SaddleSystem(LinearMap{T}(A),LinearMap{T}(B₂),LinearMap{T}(B₁ᵀ),
                        LinearMap{T}(C),Ainv,linear_map(filter,zeros(T,size(C,1)),eltype=T),
                        Vector{T},Vector{T};solver=solver,kwargs...)
end

# For cases in which C is zero, no need to pass along the argument
SaddleSystem(A::AbstractMatrix{T},B₂::AbstractMatrix{T},B₁ᵀ::AbstractMatrix{T};solver::Type{TS}=Direct,filter=I,kwargs...) where {T,TS<:SchurSolverType} =
        SaddleSystem(A,B₂,B₁ᵀ,zeros(T,size(B₂,1),size(B₁ᵀ,2));solver=solver,filter=filter,kwargs...)

### Operators are functions or function-like operators
# This version should take in functions or function-like objects that act upon given
# data types u and f. Should transform them into operators that act on abstract vectors
# of the same size
# There should already be an \ operator associated with A
# NOTE: should change default value of eltype to eltype(u)
function SaddleSystem(A,B₂,B₁ᵀ,C,u::TU,f::TF;eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TU,TF,TS<:SchurSolverType}

    return SaddleSystem(linear_map(A,u,eltype=eltype),linear_map(B₂,u,f,eltype=eltype),
                        linear_map(B₁ᵀ,f,u,eltype=eltype),
                        linear_map(C,f,eltype=eltype),
                        linear_inverse_map(A,u,eltype=eltype),
                        linear_map(filter,f,eltype=eltype),TU,TF;solver=solver,kwargs...)
end

# No C operator provided, so set it to zero
SaddleSystem(A,B₂,B₁ᵀ,u::TU,f::TF;
            eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TU,TF,TS<:SchurSolverType} =
            SaddleSystem(A,B₂,B₁ᵀ,C_zero(f,eltype),u,f;eltype=eltype,filter=filter,solver=solver,kwargs...)

# Unconstrained system
SaddleSystem(A,u::TU;eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TU,TS<:SchurSolverType} = SaddleSystem(A,nothing,nothing,u,Type{eltype}[];eltype=eltype,filter=filter,solver=solver,kwargs...)

### For handling ArrayPartition arguments for the solution/rhs
SaddleSystem(A,B₂,B₁ᵀ,C,v::ArrayPartition;
            eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TS<:SchurSolverType} =
            SaddleSystem(A,B₂,B₁ᵀ,C,v.x[1],v.x[2];eltype=eltype,filter=filter,solver=solver,kwargs...)

# No C operator
SaddleSystem(A,B₂,B₁ᵀ,v::ArrayPartition;
            eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TS<:SchurSolverType} =
            SaddleSystem(A,B₂,B₁ᵀ,v.x[1],v.x[2];eltype=eltype,filter=filter,solver=solver,kwargs...)

# Unconstrained system
SaddleSystem(A,v::ArrayPartition;
            eltype=Float64,filter=I,solver::Type{TS}=Direct,kwargs...) where {TS<:SchurSolverType} =
            SaddleSystem(A,v.x[1];eltype=eltype,filter=filter,solver=solver,kwargs...)


### AUXILIARY ROUTINES

function Base.show(io::IO, S::SaddleSystem{T,Ns,Nc,TU,TF,TS}) where {T,Ns,Nc,TU,TF,TS<:SchurSolverType}
     println(io, "Saddle system with $Ns states and $Nc constraints and")
     println(io, "   State vector of type $TU")
     println(io, "   Constraint vector of type $TF")
     println(io, "   Elements of type $T")
     println(io, "using a $TS solver")
end

"""
    Base.size(::SaddleSystem)

Report the size of a [`SaddleSystem`](@ref).
"""
size(::SaddleSystem{T,Ns,Nc}) where {T,Ns,Nc} = (Ns+Nc,Ns+Nc)

"""
    Base.eltype(::SaddleSystem)

Report the element type of a [`SaddleSystem`](@ref).
"""
eltype(::SaddleSystem{T,Ns,Nc}) where {T,Ns,Nc} = T

C_zero(f,eltype) = zeros(eltype,length(f),length(f))

function _check_sizes(A,B₂,B₁ᵀ,C,P)
    mA, nA = size(A)
    mB1, nB1 = size(B₁ᵀ)
    mB2, nB2 = size(B₂)
    mC, nC = size(C)
    mP, nP = size(P)

    # check compatibility of sizes
    mA == nA  || error("A is not square")
    mA == mB1 || error("Incompatible number of rows in A and B₁ᵀ")
    nA == nB2 || error("Incompatible number of columns in A and B₂")
    mC == mB2 || error("Incompatible number of rows in C and B₂")
    nC == nB1 || error("Incompatible number of columns in C and B₁ᵀ")
    mP == nP == mC  || error("Filter has incompatible dimensions")

    ns = nA
    nc = nB1

    return ns, nc
end
