include("../src/GeometricAlgebra.jl")

project(a::MVector, B::MVector) = dot(a, B) * inv(B)
reject(a::MVector, B::MVector) = ∧(a, B) * inv(B)

reflect(a::MVector, U::MVector) = (u = conj(U); -inv(normsq(u)) * u * a * inv(u))