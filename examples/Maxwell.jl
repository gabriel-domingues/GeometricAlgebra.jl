include("../src/GeometricAlgebra.jl")

using Symbolics
using LinearAlgebra

@variables Ex Ey Ez Bx By Bz ρ Jx Jy Jz c μ0 ϵ0 ∂t ∂x ∂y ∂z

μ0 = 1 / (c^2 * ϵ0)

MaxwellEq(s::Symbol) = MaxwellEq(Val(s))

## Vector calculus

divg(s::Vector) = dot([∂x, ∂y, ∂z], s)
curl(s::Vector) = cross([∂x, ∂y, ∂z], s)

function MaxwellEq(::Val{:vector})
    E = [Ex, Ey, Ez]
    B = [Bx, By, Bz]
    J = [Jx, Jy, Jz]

    GaussELaw = divg(E) .- ρ / ϵ0
    AmpereLaw = (∂t / c) .* E - c .* curl(B) .+ J .* inv(c * ϵ0)
    FaradayLaw = curl(E) .+ ∂t .* B
    GaussMLaw = c .* divg(B)
    [
        GaussELaw,
        AmpereLaw[1],
        AmpereLaw[2],
        FaradayLaw[3],
        AmpereLaw[3],
        -FaradayLaw[2],
        FaradayLaw[1],
        GaussMLaw
    ] .~ 0 .|> simplify
end

## Geometric Algebra

function MaxwellEq(::Val{:geometric})
    i = e1 * e2 * e3

    ∇ = (∂t / c) * 1 + ∂x * e1 + ∂y * e2 + ∂z * e3
    E = Ex * e1 + Ey * e2 + Ez * e3
    B = Bx * e1 + By * e2 + Bz * e3
    F = E + i * c * B
    J = (c * ρ) * 1 - Jx * e1 - Jy * e2 - Jz * e3
    (∇ * F - J * inv(c * ϵ0)).t .~ 0 .|> simplify
end

println("Vector equations:")
display(MaxwellEq(:vector))
println("Geometric equations:")
display(MaxwellEq(:geometric))