using Interpolations

export interpolatable_field

"""
    interpolatable_field(x,y,f::ScalarGridData)

Generates an interpolatable version of grid data `f`, based on coordinates in
`x` and `y` (which should be in range form). The output can be called as a function
with coordinate pairs as arguments.
"""
function interpolatable_field(x::AbstractArray,y::AbstractArray,f::ScalarGridData)

  (length(x) == size(f,1) && length(y) == size(f,2)) ||
    error("incompatible sizes of field and coordinate data")

  return CubicSplineInterpolation((x, y),f, extrapolation_bc = (Flat(),Flat()))

end

"""
    interpolatable_field(xu,yu,xv,yv,f::EdgeData/NodePair)

Generates an interpolatable version of grid data `f`, based on coordinates in
`xu`, `yu`, `xv`, `yv` (which should be in range form). The output can be called as a function
with coordinate pairs as arguments.
"""
function interpolatable_field(xu::AbstractRange,yu::AbstractRange,
                              xv::AbstractRange,yv::AbstractRange,q::Union{Edges,NodePair})

  return interpolatable_field(xu,yu,q.u), interpolatable_field(xv,yv,q.v)

end

"""
    interpolatable_field(f::GridData,g::PhysicalGrid)

Generates an interpolatable version of grid data `f`, based on grid `g`. The output can be called as a function
with coordinate pairs as arguments.
"""
interpolatable_field(f::GridData,g::PhysicalGrid) = interpolatable_field(coordinates(f,g)...,f)
