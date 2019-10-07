module RigidBodyMotions

export RigidBodyMotion, Kinematics, d_dt, assign_velocity!

using DocStringExtensions
import ForwardDiff
import Base: +, *, -, >>, <<, show

using Compat
using Compat: round

"""
An abstract type for types that takes in time and returns `(c, ċ, c̈, α, α̇, α̈)`.
"""
abstract type Kinematics end

"""
An abstract type for real-valued functions of time.
"""
abstract type Profile end

"""
    RigidBodyMotion

A type to store the body's current kinematics

# Fields

- `c`: current centroid position (relative to initial position)
- `ċ`: current centroid velocity
- `c̈`: current centroid acceleration
- `α`: current angle (relative to initial angle)
- `α̇`: current angular velocity
- `α̈`: current angular acceleration
- `kin`: a [`Kinematics`](@ref) structure

The first six fields are meant as a cache of the current kinematics
while the `kin` field can be used to find the plate kinematics at any time.
"""
mutable struct RigidBodyMotion
    c::ComplexF64
    ċ::ComplexF64
    c̈::ComplexF64
    α::Float64
    α̇::Float64
    α̈::Float64

    kin::Kinematics
end

RigidBodyMotion(ċ, α̇) = RigidBodyMotion(0.0im, complex(ċ), 0.0im, 0.0, float(α̇),
                                          0.0, Constant(ċ, α̇))
RigidBodyMotion(kin::Kinematics) = RigidBodyMotion(kin(0)..., kin)
(m::RigidBodyMotion)(t) = m.kin(t)


function (m::RigidBodyMotion)(t,x̃::Tuple{Real,Real})
  # This expects coordinates in body's own coordinate system
  #
  z̃ = ComplexF64(x̃[1],x̃[2])
  m.c, m.ċ, m.c̈, m.α, m.α̇, m.α̈ = m.kin(t)
  z = exp(im*m.α)*z̃
  return m.c + z, m.ċ + im*m.α̇*z, m.c̈ + (im*m.α̈-m.α̇^2)*z
end

"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     x::AbstractVector{Float64},y::AbstractVector{Float64},
                     xc::Real,yc::Real,α::Real,
                     motion,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
function assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                          x::AbstractVector{Float64},y::AbstractVector{Float64},
                          xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real)
   _,ċ,_,_,α̇,_ = m(t)
  for i = 1:length(x)
      Δz = (x[i]-xc)+im*(y[i]-yc)
      ċi = ċ + im*α̇*Δz
      u[i] = real(ċi)
      v[i] = imag(ċi)
  end
  nothing
end


function show(io::IO, m::RigidBodyMotion)
    println(io, "Rigid Body Motion:")
    println(io, "  ċ = $(round(m.ċ, digits=2))")
    println(io, "  c̈ = $(round(m.c̈, digits=2))")
    println(io, "  α̇ = $(round(m.α̇, digits=2))")
    println(io, "  α̈ = $(round(m.α̈, digits=2))")
    print(io, "  $(m.kin)")
end

#=
Kinematics
=#


struct Constant{C <: Complex, A <: Real} <: Kinematics
    ċ::C
    α̇::A
end
Constant(ċ, α̇) = Constant(complex(ċ), α̇)
(c::Constant{C})(t) where C = zero(C), c.ċ, zero(C), 0.0, c.α̇, 0.0
show(io::IO, c::Constant) = print(io, "Constant (ċ = $(c.ċ), α̇ = $(c.α̇))")

"""
    Pitchup <: Kinematics

Kinematics describing a pitchup motion (horizontal translation with rotation)

# Constructors
# Fields
$(FIELDS)
"""
struct Pitchup <: Kinematics
    "Freestream velocity"
    U₀::Float64
    "Axis of rotation, relative to the plate centroid"
    a::Float64

    "Non-dimensional pitch rate ``K = \\dot{\\alpha}_0\\frac{c}{2U_0}``"
    K::Float64

    "Initial angle of attack"
    α₀::Float64
    "Nominal start of pitch up"
    t₀::Float64

    "Total pitching angle"
    Δα::Float64

    α::Profile
    α̇::Profile
    α̈::Profile
end

function Pitchup(U₀, a, K, α₀, t₀, Δα, ramp)
    Δt = 0.5Δα/K
    p = ConstantProfile(α₀) + 2K*((ramp >> t₀) - (ramp >> (t₀ + Δt)))
    ṗ = d_dt(p)
    p̈ = d_dt(ṗ)
    Pitchup(U₀, a, K, α₀, t₀, Δα, c, ċ, c̈, p, ṗ, p̈)
end

function (p::Pitchup)(t)
    α = p.α(t)
    α̇ = p.α̇(t)
    α̈ = p.α̈(t)

    c = p.U₀*t - p.a*exp(im*α)
    ċ = p.U₀ - p.a*im*α̇*exp(im*α)
    if (t - p.t₀) > p.Δα/p.K
        c̈ = 0.0im
    else
        c̈ = p.a*exp(im*α)*(α̇^2 - im*α̈)
    end

    return c, ċ, c̈, α, α̇, α̈
end

function show(io::IO, p::Pitchup)
    print(io, "Pitch-up kinematics with rate K = $(p.K)")
end


"""
    PitchHeave <: Kinematics

Kinematics describing an oscillatory pitching and heaving (i.e. plunging) motion

# Constructors
# Fields
$(FIELDS)
"""
struct PitchHeave <: Kinematics
    "Freestream velocity"
    U₀::Float64

    "Axis of pitch rotation, relative to the plate centroid"
    a::Float64

    "Reduced frequency ``K = \\frac{\\Omega c}{2U_0}``"
    K::Float64

    "Phase of pitch (in radians)"
    ϕp::Float64

    "Phase of heave (in radians)"
    ϕh::Float64

    "Mean angle of attack"
    α₀::Float64

    "Amplitude of pitching"
    Δα::Float64

    "Amplitude of translational heaving"
    A::Float64

    Y::Profile
    Ẏ::Profile
    Ÿ::Profile

    α::Profile
    α̇::Profile
    α̈::Profile
end

function PitchHeave(U₀, a, K, ϕp, α₀, Δα, A, ϕh)
    p = A*(Sinusoid(2K) >> (ϕp/(2K)))
    ṗ = d_dt(p)
    p̈ = d_dt(ṗ)
    α = ConstantProfile(α₀) + Δα*(Sinusoid(2K) >> (ϕh/(2K)))
    α̇ = d_dt(α)
    α̈ = d_dt(α̇)
    PitchHeave(U₀, a, K, ϕp, ϕh, α₀, Δα, A, p, ṗ, p̈, α, α̇, α̈)
end

function (p::PitchHeave)(t)
    α = p.α(t)
    α̇ = p.α̇(t)
    α̈ = p.α̈(t)

    c = p.U₀*t + im*p.Y(t) - p.a*exp(im*α)
    ċ = p.U₀ + im*p.Ẏ(t) - p.a*im*α̇*exp(im*α)
    c̈ = im*p.Ÿ(t) + p.a*exp(im*α)*(α̇^2 - im*α̈)

    return c, ċ, c̈, α, α̇, α̈
end

function show(io::IO, p::PitchHeave)
    println(io, "Oscillatory pitch-heave kinematics with")
    println(io, "     Reduced frequency K = $(p.K)")
    println(io, "     Heaving amplitude A = $(p.A)")
    println(io, "     Pitching amplitude Δα = $(p.Δα)")
    println(io, "     Pitch lag ϕp = $(p.ϕp)")
    println(io, "     Heave lag ϕh = $(p.ϕh)")
end

struct Oscillation <: Kinematics
    "Angular frequency"
    Ω :: Float64

    "Amplitude x direction"
    Ax:: Float64

    "Phase in x direction."
    ϕx :: Float64

    "Amplitude y direction"
    Ay:: Float64

    "Phase in y direction."
    ϕy :: Float64

    cx::Profile
    ċx::Profile
    c̈x::Profile

    cy::Profile
    ċy::Profile
    c̈y::Profile

end

function Oscillation(Ω,Ax,ϕx,Ay,ϕy)
    Δtx = ϕx/Ω
    px = Ax*(Sinusoid(Ω) << Δtx)
    ṗx = d_dt(px)
    p̈x = d_dt(ṗx)

    Δty = ϕy/Ω
    py = Ay*(Sinusoid(Ω) << Δty)
    ṗy = d_dt(py)
    p̈y = d_dt(ṗy)
    Oscillation(Ω, Ax, ϕx, Ay, ϕy, px, ṗx, p̈x, py, ṗy, p̈y)
end

function (p::Oscillation)(t)
    α = 0.0
    α̇ = 0.0
    α̈ = 0.0

    c = ComplexF64(p.cx(t)) + im*ComplexF64(p.cy(t))
    ċ = ComplexF64(p.ċx(t)) + im*ComplexF64(p.ċy(t))
    c̈ = ComplexF64(p.c̈x(t)) + im*ComplexF64(p.c̈y(t))
    return c, ċ, c̈, α, α̇, α̈

    #return [p.ċ(t),0.0], [p.c̈(t),0.0], α̇
end

struct OscilX <: Kinematics
    "Angular frequency"
    Ω :: Float64

    "Mean velocity"
    Umean :: Float64

    "Velocity amplitude"
    Ux:: Float64

    "Velocity phase"
    ϕx :: Float64

    cx::Profile
    ċx::Profile
    c̈x::Profile

end

function OscilX(Ω,Umean,Ux,ϕx)
    Δtx = ϕx/Ω
    px = ConstantProfile(0.0)
    ṗx = ConstantProfile(Umean) + Ux*(Sinusoid(Ω) << Δtx)
    p̈x = d_dt(ṗx)

    OscilX(Ω, Umean, Ux, ϕx, px, ṗx, p̈x)
end

function (p::OscilX)(t)
    α = 0.0
    α̇ = 0.0
    α̈ = 0.0

    c = ComplexF64(p.cx(t))
    ċ = ComplexF64(p.ċx(t))
    c̈ = ComplexF64(p.c̈x(t))
    return c, ċ, c̈, α, α̇, α̈

    #return [p.ċ(t),0.0], [p.c̈(t),0.0], α̇
end

abstract type Switch end
abstract type SwitchOn <: Switch end
abstract type SwitchOff <: Switch end

"""
    SwitchedKinematics <: Kinematics

Modulates a given set of kinematics between simple on/off states. The velocity
specified by the given kinematics is toggled on/off.

# Fields
$(FIELDS)
"""
struct SwitchedKinematics{S <: Switch} <: Kinematics

    "time at which the kinematics should be turned on"
    t_on :: Float64

    "time at which the kinematics should be turned off"
    t_off :: Float64

    "kinematics to be followed in the on state"
    kin :: Kinematics

    off :: Kinematics

    SwitchedKinematics(t_on,t_off,kin) = t_on > t_off ?
            new{SwitchOn}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0)) :
            new{SwitchOff}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0))
end

# note that these do not introduce impulsive changes into the derivatives
(p::SwitchedKinematics{SwitchOn})(t) = t <= p.t_on ? p.off(t) : p.kin(t-p.t_on)

(p::SwitchedKinematics{SwitchOff})(t) = t <= p.t_off ? p.kin(t-p.t_on) : p.off(t)


#=
Profiles
=#


"""
    ConstantProfile(c::Number)

Create a profile consisting of a constant `c`.

# Example

```jldoctest
julia> p = RigidBodyMotions.ConstantProfile(1.0)
Constant (2.3)
```
"""
struct ConstantProfile <: Profile
    c::Number
end

function show(io::IO, p::ConstantProfile)
    print(io, "Constant ($(p.c))")
end

(p::ConstantProfile)(t) = p.c

struct DerivativeProfile{P} <: Profile
    p::P
end

function show(io::IO, ṗ::DerivativeProfile)
    print(io, "d/dt ($(ṗ.p))")
end

(ṗ::DerivativeProfile)(t) = ForwardDiff.derivative(ṗ.p, t)

"""
    d_dt(p::Profile)

Take the time derivative of `p` and return it as a new profile.

# Example

```jldoctest
julia> s = Plates.RigidBodyMotions.Sinusoid(π)
Sinusoid (ω = 3.14)

julia> s.([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 0.0
 1.0
 0.707107

julia> c = Plates.RigidBodyMotions.d_dt(s)
d/dt (Sinusoid (ω = 3.14))

julia> c.([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
  3.14159
  1.92367e-16
 -2.22144
```
"""
d_dt(p::Profile) = DerivativeProfile(p)

struct ScaledProfile{N <: Real, P <: Profile} <: Profile
    s::N
    p::P
end
function show(io::IO, p::ScaledProfile)
    print(io, "$(p.s) × ($(p.p))")
end

"""
    s::Number * p::Profile

Returns a scaled profile with `(s*p)(t) = s*p(t)`

# Example

```jldoctest
julia> s = Plates.RigidBodyMotions.Sinusoid(π)
Sinusoid (ω = 3.14)

julia> 2s
2 × (Sinusoid (ω = 3.14))

julia> (2s).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 0.0
 2.0
 1.41421
```
"""
s::Number * p::Profile = ScaledProfile(s, p)

"""
    -(p₁::Profile, p₂::Profile)

```jldoctest
julia> s = Plates.RigidBodyMotions.Sinusoid(π)
Sinusoid (ω = 3.14)

julia> 2s
2 × (Sinusoid (ω = 3.14))

julia> (2s).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 0.0
 2.0
 1.41421

julia> s = Plates.RigidBodyMotions.Sinusoid(π);

julia> s.([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 0.0
 1.0
 0.707107

julia> (-s).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 -0.0
 -1.0
 -0.707107

julia> (s - s).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
```
"""
-(p::Profile) = ScaledProfile(-1, p)
(p::ScaledProfile)(t) = p.s*p.p(t)

struct ShiftedProfile{N <: Real, P <: Profile} <: Profile
    Δt::N
    p::P
end
function show(io::IO, p::ShiftedProfile)
    print(io, "$(p.p) >> $(p.Δt)")
end

(p::ShiftedProfile)(t) = p.p(t - p.Δt)

"""
    p::Profile >> Δt::Number

Shift the profile in time so that `(p >> Δt)(t) = p(t - Δt)`

# Example

```jldoctest
julia> s = Plates.RigidBodyMotions.Sinusoid(π);

julia> s >> 0.5
Sinusoid (ω = 3.14) >> 0.5

julia> (s >> 0.5).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
 -1.0
  0.0
  0.707107

julia> (s << 0.5).([0.0, 0.5, 0.75])
3-element Array{Float64,1}:
  1.0
  1.22465e-16
 -0.707107
```
"""
p::Profile >> Δt::Number = ShiftedProfile(Δt, p)
p::Profile << Δt::Number = ShiftedProfile(-Δt, p)

struct AddedProfiles{T <: Tuple} <: Profile
    ps::T
end
function show(io::IO, Σp::AddedProfiles)
    println(io, "AddedProfiles:")
    for p in Σp.ps
        println(io, "  $p")
    end
end

"""
    p₁::Profile + p₂::Profile

Add the profiles so that `(p₁ + p₂)(t) = p₁(t) + p₂(t)`.

# Examples

```jldoctest
julia> ramp₁ = Plates.RigidBodyMotions.EldredgeRamp(5)
logcosh ramp (aₛ = 5.0)

julia> ramp₂ = Plates.RigidBodyMotions.ColoniusRamp(5)
power series ramp (n = 5.0)

julia> ramp₁ + ramp₂
AddedProfiles:
  logcosh ramp (aₛ = 5.0)
  power series ramp (n = 5.0)


julia> ramp₁ + (ramp₂ + ramp₁) == ramp₁ + ramp₂ + ramp₁
true

```
"""
+(p::Profile, Σp::AddedProfiles) = AddedProfiles((p, Σp.ps...))
+(Σp::AddedProfiles, p::Profile) = AddedProfiles((Σp.ps..., p))
function +(Σp₁::AddedProfiles, Σp₂::AddedProfiles)
    AddedProfiles((Σp₁..., Σp₂...))
end

-(p₁::Profile, p₂::Profile) = p₁ + (-p₂)

+(p::Profile...) = AddedProfiles(p)

function (Σp::AddedProfiles)(t)
    f = 0.0
    for p in Σp.ps
        f += p(t)
    end
    f
end


struct MultipliedProfiles{T <: Tuple} <: Profile
    ps::T
end
function show(io::IO, Πp::MultipliedProfiles)
    println(io, "MultipliedProfiles:")
    for p in Πp.ps
        println(io, "  $p")
    end
end
*(p::Profile, Πp::MultipliedProfiles) = MultipliedProfiles((p, Πp.ps...))
*(Πp::MultipliedProfiles, p::Profile) = MultipliedProfiles((Πp.ps..., p))
function *(Πp₁::MultipliedProfiles, Πp₂::MultipliedProfiles)
    MultipliedProfiles((Πp₁..., Πp₂...))
end

function (Πp::MultipliedProfiles)(t)
    f = 1.0
    for p in Πp.ps
        f *= p(t)
    end
    f
end


struct Sinusoid <: Profile
    ω::Float64
end
(s::Sinusoid)(t) = sin(s.ω*t)
show(io::IO, s::Sinusoid) = print(io, "Sinusoid (ω = $(round(s.ω, digits=2)))")

struct EldredgeRamp <: Profile
    aₛ::Float64
end
(r::EldredgeRamp)(t) = 0.5(log(2cosh(r.aₛ*t)) + r.aₛ*t)/r.aₛ
show(io::IO, r::EldredgeRamp) = print(io, "logcosh ramp (aₛ = $(round(r.aₛ, digits=2)))")

struct ColoniusRamp <: Profile
    n::Int
end
function (r::ColoniusRamp)(t)
    Δt = t + 0.5
    if Δt ≤ 0
        0.0
    elseif Δt ≥ 1
        Δt - 0.5
    else
        f = 0.0
        for j = 0:r.n
            f += binomial(r.n + j, j)*(r.n - j + 1)*(1 - Δt)^j
        end
        f*Δt^(r.n + 2)/(2r.n + 2)
    end
end
show(io::IO, r::ColoniusRamp) = print(io, "power series ramp (n = $(round(r.n, digits=2)))")

end
