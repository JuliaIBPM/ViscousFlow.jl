import Base: *, \

# Identity operator

struct Identity end

(*)(::Identity,s::GridData) = s

# Lots of other operators

include("innerproducts.jl")
include("convolution.jl")
include("lgf.jl")
include("ddf.jl")
include("laplacian.jl")
include("intfact.jl")
include("diffcalculus.jl")
include("shift.jl")
include("regularization.jl")
