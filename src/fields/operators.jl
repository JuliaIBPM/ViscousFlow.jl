include("innerproducts.jl")
include("convolution.jl")
include("lgf.jl")
include("ddf.jl")

import Base: *, \

# Identity

struct Identity end

(*)(::Identity,s::Union{Nodes,Edges}) = s

include("laplacian.jl")
include("intfact.jl")
include("diffcalculus.jl")
include("differencing1d.jl")
include("shift.jl")
include("regularization.jl")
