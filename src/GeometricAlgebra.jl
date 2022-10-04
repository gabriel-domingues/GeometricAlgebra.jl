module GeometricAlgebra

export
    MVector,
    dim,
    signature,
    scalar,
    rev,
    MComplex,
    MQuaternion,
    STA,
    e1,
    e2,
    e3,
    γ0,
    γ1,
    γ2,
    γ3

include("GA.jl")
include("names.jl")
include("dispatch.jl")

end