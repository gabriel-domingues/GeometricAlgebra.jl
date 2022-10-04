include("../src/GeometricAlgebra.jl")

struct Lorentz{T<:Real}
    γ0_::STA{T}
    norm::T
    Lorentz(γ::STA{T}) where {T<:Real} = new{T}(γ, norm(γ))
end

Base.show(io::IO, L::Lorentz{T}) where {T} = print(io, "Lorentz(", L.γ0_::STA{T}, ")")

(L::Lorentz{T})(m::STA{T}) where {T<:Real} = m * L.γ0_ * γ0 * inv(L.norm)