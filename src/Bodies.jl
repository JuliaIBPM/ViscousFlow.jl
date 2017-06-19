module Bodies

export Body

import Whirl2d
import Whirl2d.Grids


struct BodyConfig

    "Position of body reference point in inertial space"
    xref::Array{Float64,1}

    "Rotation tensor of body"
    rot::Array{Float64,2}

end

function BodyConfig(xref::Vector{<:Real},angle::Float64)
    rot = [[cos(angle) -sin(angle)]; [sin(angle) cos(angle)]]

    BodyConfig(xref,rot)

end

transform(x::Array{Float64,1},config::BodyConfig) = config.xref + config.rot*x

function transform(x::Array{Array{Float64,1},1},config::BodyConfig)
    [transform(x[i],config) for i=1:length(x)]
end


mutable struct Body
    "number of Lagrange points on body"
    N::Int

    "coordinates of Lagrange points in body-fixed system"
    xtilde::Array{Array{Float64,1},1}

    "coordinates of Lagrange points in inertial system"
    x::Array{Array{Float64,1},1}

    "Lagrange forcing field components"
    f::Array{Array{Float64,1},1}

    "body velocity components at Lagrange points"
    vel::Array{Array{Float64,1},1}

    "body configuration"
    config::BodyConfig

end

# Set up blank body
Body() = Body(0,[])

function Body(N,xtilde)

    # set up array of inertial coordinates
    x = xtilde

    # set up array of forcing field
    f = [zeros(Whirl2d.ndim) for i=1:N]

    # set up zero array of body velocities
    vel = [zeros(Whirl2d.ndim) for i=1:N]

    # set configuration to origin and zero angle
    config = BodyConfig([0.0,0.0],0.0);

    Body(N,xtilde,x,f,vel,config)

end

function Body(N,xtilde,config)

    b = Body(N,xtilde)
    update_body!(b,config)
    b

end

Body(N,xtilde,xref,angle) = Body(N,xtilde,BodyConfig(xref,angle))

function dims(body::Body)
    xmin = Inf*ones(Whirl2d.ndim)
    xmax = -Inf*ones(Whirl2d.ndim)
    for x in body.x
        xmin = [min(xmin[j],x[j]) for j = 1:Whirl2d.ndim]
        xmax = [max(xmax[j],x[j]) for j = 1:Whirl2d.ndim]
    end
    xmin, xmax
end

function update_body!(body::Body,config::BodyConfig)
    body.config = config
    body.x = transform(body.xtilde,config)
end

function set_body_velocities!(body::Body,vel::Array{Array{Float64,1},1})
    body.vel = vel
end

function Circle(N::Int,rad)::Body

    # set up the points on the circle with radius `rad`
    x = [[rad*cos(2*pi*(i-1)/N),rad*sin(2*pi*(i-1)/N)] for i=1:N]

    # put it at the origin, with zero angle
    Body(N,x,[0.0,0.0],0.0)

end

function Circle(N::Int,rad,xcent::Vector{<:Real},angle)::Body

    b = Circle(N,rad)
    update_body!(b,BodyConfig(xcent,angle))
    b

end

function Plate(N::Int,len)::Body

    # set up points on plate
    x = [[len*(-0.5 + 0.5*(i-1)/(N-1)),0.0] for i=1:N]

    # put it at the origin, with zero angle
    Body(N,x,[0.0,0.0],0.0)

end

function Plate(N::Int,rad,xcent::Vector{<:Real},angle)::Body

    b = Plate(N,rad)
    update_body!(b,BodyConfig(xcent,angle))
    b

end


function Base.show(io::IO, b::Body)
    println(io, "Body: number of points = $(b.N), "*
    		"reference point = ($(b.config.xref[1]),$(b.config.xref[2])), "*
		"rotation matrix = $(b.config.rot)")
end


end
