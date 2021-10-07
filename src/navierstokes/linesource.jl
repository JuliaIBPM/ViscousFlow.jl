struct LineSourceParams{BT<:Body}
    body :: BT
end

function LineSourceParams(x::Vector{T},y::Vector{T};closure=:closed) where {T}
      closuretype = (closure == :closed ? RigidBodyTools.ClosedBody : RigidBodyTools.OpenBody)
      LineSourceParams(BasicBody(x,y,closuretype=closuretype))
end

struct PrescribedLineSource{BT<:Body,AT,GT<:ScalarGridData,ST<:ScalarData,RT<:RegularizationMatrix,ET<:InterpolationMatrix}
    body :: BT
    arccoord :: AT
    cache1 :: GT
    q :: ST
    R :: RT
    E :: ET
end

function PrescribedLineSource(params::LineSourceParams,u::ScalarGridData,g::PhysicalGrid;ddftype=CartesianGrids.Yang3)
    @unpack body = params
    pts = VectorData(collect(body))
    regop = _regularization(pts,g,body,ddftype)

    s = ScalarData(pts)

    arccoord = ScalarData(pts)
    arccoord .= accumulate(+,dlengthmid(body))

    R = RegularizationMatrix(regop,s,u)
    E = InterpolationMatrix(regop,u,s)
    PrescribedLineSource(body,arccoord,similar(u),similar(s),R,E)
end

set_linesource_strength!(qline::PrescribedLineSource{BT,GT,ST},q::ST) where {BT,GT,ST} = qline.q .= q

function set_linesource_strength!(qline::PrescribedLineSource,q::AbstractVector{T}) where {T<:Real}
    @assert length(q) == length(qline.q)
    set_linesource_strength!(qline,ScalarData(q))
end
