struct LineSourceParams{BT<:Body}
    body :: BT
end

function LineSourceParams(x::Vector{T},y::Vector{T};closure=:closed) where {T}
      closuretype = (closure == :closed ? RigidBodyTools.ClosedBody : RigidBodyTools.OpenBody)
      LineSourceParams(BasicBody(x,y,closuretype=closuretype))
end

struct PrescribedLineSource{BT<:Body,AT,GT<:VectorGridData,ST<:VectorData,RT<:RegularizationMatrix,ET<:InterpolationMatrix}
    body :: BT
    arccoord :: AT
    cache1 :: GT
    τ :: ST
    R :: RT
    E :: ET
end

function PrescribedLineSource(params::LineSourceParams,u::VectorGridData,g::PhysicalGrid;ddftype=CartesianGrids.Yang3)
    @unpack body = params
    
    pts = VectorData(collect(body))
    regop = _regularization(pts,g,body,ddftype)

    s = VectorData(pts)

    arccoord = ScalarData(pts)
    arccoord .= accumulate(+,dlengthmid(body))

    R = RegularizationMatrix(regop,s,u)
    E = InterpolationMatrix(regop,u,s)
    PrescribedLineSource(body,arccoord,similar(u),similar(s),R,E)
end

set_linesource_strength!(τ_bc::PrescribedLineSource{BT,AT,GT,ST},q::ST) where {BT,AT,GT,ST} = τ_bc.τ .= q

# function set_linesource_strength!(τ_bc::PrescribedLineSource,q::AbstractVector{T}) where {T<:Real}
#     @assert length(q) == length(τ_bc.τ)
#     set_linesource_strength!(τ_bc,VectorData(q))
# end