### Right-hand side and solution vectors

const SaddleVector = ArrayPartition

#SaddleVector(u::TU,f::TF) where {TU,TF} = ArrayPartition(u,f)

state(u::SaddleVector) = u.x[1]
constraint(u::SaddleVector) = u.x[2]
