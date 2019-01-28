
struct PointForce{T}
    x :: Tuple{Float64,Float64}
    f0 :: Union{Float64,Tuple{Float64,Float64}}
    t0 :: Float64
    σ :: Float64
    fbuffer :: Union{ScalarData,VectorData}
    ubuffer :: T
    regop :: Regularize
end

"""
    PointForce(u::Union{Nodes,Edges},x0::Tuple{Float64,Float64},f0,t0,σ,sys::NavierStokes)

Constructor function that immerses a point force in the `u`-type data of system `sys`,
of strength `f0` to be applied at physical position `x0`, modulated by a Gaussian
centered at time `t0` with standard deviation `σ`. The data `u` should be of either
`Nodes` or `Edges` type. If `Nodes`, then `f0` should be a scalar; if `Edges`,
then `f0` should be a tuple.

The resulting function is a function of time and generates a field on `u`-type data.
"""
function PointForce(w₀::T,x0::Tuple{Float64,Float64},
                    f0::Union{Float64,Tuple{Float64,Float64}},
                    t0::Float64,σ::Float64,
        sys::NavierStokes) where {T <: Union{Nodes,Edges}}
    Xf = VectorData([x0[1]],[x0[2]])
    Ff = (T <: Nodes) ? ScalarData(Xf) : VectorData(Xf)
    Ff .= f0

    regop = Regularize(Xf,Fields.cellsize(sys);I0=Fields.origin(sys),issymmetric=true)
    PointForce{T}(x0,f0,t0,σ,Ff,T(),regop)
end

(f::PointForce)(t) = f.regop(f.ubuffer,rmul!(deepcopy(f.fbuffer),exp(-(t-f.t0)^2/f.σ^2)))

function Base.show(io::IO,f::PointForce{T}) where {T}
    println(io,"Transient point force applied on the $(T) field.")
    println(io,"   strength = $(f.f0)")
    println(io,"   location = $(f.x)")
    println(io,"   central time = $(f.t0)")
    println(io,"   half-interval = $(f.σ)")
end
