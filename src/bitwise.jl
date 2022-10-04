function round2(x::Integer)
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x |= x >> 32
    x += 1
end

@inline bitlog2(x::Integer) = count_ones(round2(x) - 1)
@inline parity(x::Integer) = 1 - ((x & 1) << 1)

@inline bin2graycode(x::Integer) = x ⊻ (x >> 1)
