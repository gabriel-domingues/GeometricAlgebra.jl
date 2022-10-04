include("bitwise.jl")

extend(t::Vector{T}, n::Int64) where {T} = vcat(t, zeros(T, max(n - length(t), 0)))

struct MVector{Sign,T<:Real} <: Number
    t::Vector{T}

    MVector{T}() where {T<:Real} = new{(0, 0),T}([])

    MVector{Sign}(t::Vector{T}) where {Sign,T<:Real} = new{Sign,T}(extend(t, round2(length(t))))
    MVector{Sign}(t::Vararg{T}) where {Sign,T<:Real} = MVector{Sign}(collect(t))

    MVector(t::Vector{T}) where {T<:Real} = (N = bitlog2(length(t)); new{(N, 0),T}(extend(t, round2(length(t)))))
    MVector(t::Vararg{T}) where {T<:Real} = MVector(collect(t))
end

dim(m::MVector) = length(m.t)
signature(::MVector{Sign}) where {Sign} = Sign

function match(a::MVector, b::MVector)
    ℓ = max(dim(a), dim(b))
    S = max(signature(a), signature(b))
    MVector{S}(extend(a.t, ℓ)), MVector{S}(extend(b.t, ℓ))
end

basis_vector(::Type{T}, n::Integer) where {T<:Real} = (a = zeros(T, round2(n + 1)); a[n+1] = one(T); a)
basis(::Type{T}, n::Integer) where {T<:Real} = MVector(basis_vector(T, n))
basis(Sign::Tuple{Int,Int}, ::Type{T}, n::Integer) where {T<:Real} = MVector{Sign}(basis_vector(T, n))

Base.show(io::IO, m::MVector{Sign,T}) where {Sign,T} = print(io, "MVector", Sign, m.t::Vector{T})

Base.zero(a::MVector{Sign,T}) where {Sign,T<:Real} = MVector{Sign}(zeros(T, dim(a)))
pseudo(a::MVector{Sign,T}) where {Sign,T<:Real} = (I = zero(a); I.t[end] = one(T); I)

Base.convert(::Type{MVector{Sign,T}}, m::MVector) where {Sign,T<:Real} = MVector{Sign}(convert(Vector{T}, m.t))
Base.convert(::Type{MVector{Sign,T}}, v::T) where {Sign,T<:Real} = MVector{Sign}(v)

Base.promote_rule(::Type{MVector{Sign,T}}, ::Type{S}) where {Sign,T<:Real,S<:Real} = MVector{Sign,promote_type(T, S)}
Base.promote_rule(::Type{MVector{Sign,T}}, ::Type{T}) where {Sign,T<:Real} = MVector{Sign,T}
Base.promote_rule(::Type{MVector{Sign1,T}}, ::Type{MVector{Sign2,S}}) where {Sign1,Sign2,T<:Real,S<:Real} =
    MVector{max(Sign1, Sign2),promote_type(T, S)}

Base.:+(a::MVector, b::MVector) = ((a, b) = match(a, b); MVector{signature(a)}(a.t .+ b.t))

Base.:-(a::MVector{Sign}) where {Sign} = MVector{Sign}(.-a.t)
Base.:-(a::MVector, b::MVector) = ((a, b) = match(a, b); MVector{signature(a)}(a.t .- b.t))

Base.:*(λ::T, m::MVector{Sign}) where {Sign,T<:Real} = MVector{Sign}(λ .* m.t)
Base.:*(m::MVector{Sign}, λ::T) where {Sign,T<:Real} = MVector{Sign}(m.t .* λ)


@inline signof(j::Int64, k::Int64) = parity(count_ones(bin2graycode(j >> 1) & k))
@inline signof(Sign, j::Int64, k::Int64) = parity(count_ones((j & k) >> Sign[1]))

function Base.:*(a::MVector, b::MVector)
    (a, b) = match(a, b)
    Sign = signature(a)
    r = 0:dim(a)-1
    c = collect(
        sum(a.t[j+1] * b.t[k+1] * signof(j, k) * signof(Sign, j, k)
            for j = r for k = r if (j ⊻ k) == i)
        for i = r)
    MVector{Sign}(c)
end

Base.:/(m::MVector{Sign}, λ::T) where {Sign,T<:Real} = MVector{Sign}(m.t ./ λ)

scalar(a::MVector) = (@assert all(iszero.(a.t[2:end])); a.t[1])

rev(a::MVector) = MVector{signature(a)}([(1 - count_ones(i - 1) & 2) * t for (i, t) in enumerate(a.t)])

Base.inv(a::MVector) = (a_ = rev(a); a_ * inv(a_ * a |> scalar))

Base.:/(a::MVector, b::MVector) = a * inv(b)
