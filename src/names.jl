const MComplex = MVector{(0, 1)}
const MQuaternion = MVector{(0, 2)}
const STA = MVector{(1, 3)}

const e1 = basis(Int, 0b001)
const e2 = basis(Int, 0b010)
const e3 = basis(Int, 0b100)

const γ0 = basis((1, 3), Int, 0b0001)
const γ1 = basis((1, 3), Int, 0b0010)
const γ2 = basis((1, 3), Int, 0b0100)
const γ3 = basis((1, 3), Int, 0b1000)
