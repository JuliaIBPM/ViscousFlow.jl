"""
    SaddlePointSystem
"""
abstract type SaddlePointSystem{T} end

struct NavierStokes <: SaddlePointSystem{false}
    domain
    A⁻¹
    L⁻¹
end

struct NavierStokesWithBody <: SaddlePointSystem{true}
    domain
    A⁻¹
    L⁻¹
    B₁ᵀ
    B₂!
    S⁻¹
    S₀⁻¹
    r₁
    r₂
end

struct TwoLevelBody <: SaddlePointSystem{true}
    domain
    A⁻¹
    L⁻¹
    B₁ᵀ
    B₂!
    S⁻¹
    S₀⁻¹
    r₁
    r₂
end
