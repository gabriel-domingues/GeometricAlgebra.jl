using Symbolics
using LinearAlgebra

Base.:*(λ::Num, m::MVector{Sign}) where {Sign} = MVector{Sign}(λ .* m.t)
Base.:*(m::MVector{Sign}, λ::Num) where {Sign} = MVector{Sign}(m.t .* λ)
Base.isless(a::MVector, b::MVector) = isless(a.t, b.t)


∧(a::MVector, b::MVector) = (a * b - rev(b) * a) / 2
LinearAlgebra.dot(a::MVector, b::MVector) = (a * b + rev(b) * a) / 2

normsq(a::MVector) = a^2 |> scalar
LinearAlgebra.norm(a::MVector) = a |> normsq |> sqrt


Base.conj(a::MVector) = a * pseudo(a)