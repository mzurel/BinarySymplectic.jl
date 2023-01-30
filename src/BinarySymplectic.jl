"""
    BinarySymplectic.jl

Tools for working with symplectic vector spaces and symplectic groups over ℤ₂.

https://github.com/mzurel/BinarySymplectic.jl
"""
module BinarySymplectic

import Base: (==), (+), (-), (*), (/), (^), (<), (>), (≤), (≥), bitstring, convert, eltype,
    hash, isequal, isless, iszero, promote, show, rand, zero
import Random: AbstractRNG, SamplerType

export AbstractSymplecticVector, SymplecticVector, SymplecticMap, Transvection, Subspace
export halfdimension, dimension, data, vector, bitstring, iszero, isequal, isless,
    symplecticform, (⋆), dotproduct, (⋅), lift, symplecticgramschmidt, transvection,
    findtransvection, symplecticgrouporder, isisotropic, isLagrangian, SYMPLECTICImproved

include("utils.jl")


############################################################################################
## Types for symplectic vectors                                                           ##
############################################################################################

abstract type AbstractSymplecticVector end

"""
    SymplecticVector{n, T}

A type that stores an element of the symplectic vector space ``ℤ₂ⁿ×ℤ₂ⁿ``.

The vector is stored as the bits in the binary expansions of a pair of integers of type `T`.
An pair of integers ``(a,b)`` represents ``a₁e₁+b₁f₁+a₂e₂+...`` where ``e₁,f₁,e₂,f₂,...``
are the standard symplectic basis vectors and ``a₁`` is the least significant bit in the
binary expansion of `a`, ``a₂`` is the second least significant bit, and so on. Similarly,
``b₁,b₂,...`` are the bits in the binary expansion of `b`.

# Examples
```jldoctest
julia> u = SymplecticVector{3, UInt8}(7, 4)
101011
```
"""
struct SymplecticVector{n, T} <: AbstractSymplecticVector
    a::T
    b::T
    function SymplecticVector{n, T}(a::T, b::T) where {n, T<:Integer}
        return new(a, b)
    end
end

"""
    typerequired(n::Integer)

Compute the smallest UInt type required to represent an element of ``ℤ₂²ⁿ``.

Elements of the symplectic vector space ``ℤ₂ⁿ×ℤ₂ⁿ`` are represented as the bits in the
binary expansions of a pair of integers. This function computes the smallest unsigned
integer type with enough bits to store an element of ``ℤ₂ⁿ``. If ``n>128``, built-in
unsigned types are too small so BigInt is returned.
"""
function typerequired(n::Integer)
    if n ≤ 128
        bitsrequired = max(8, convert(Integer, 2^ceil(log2(n))))
        return eval(Meta.parse("UInt"*string(bitsrequired)))
    else
        return BigInt
    end
end
# TODO: see if larger fixed width unsigned types from BitIntegers.jl will work for ``n>128``

"""
    SymplecticVector{n}

If the type `T` is not specified then the smallest unsigned type with enough bits to store
the vector is used.

See also [`typerequired`](@ref).

# Examples
```jldoctest
julia> v = SymplecticVector{3}(7, 4)
101011
```
"""
@inline function SymplecticVector{n}(a::Integer, b::Integer) where {n}
    T = typerequired(n)
    return SymplecticVector{n, T}(a, b)
end

# type casting
function SymplecticVector{n, T}(a::T1, b::T2) where {n, T<:Integer, T1<:Integer, T2<:Integer}
    return SymplecticVector{n, T}(convert(T, a), convert(T, b))
end

function SymplecticVector{n, T1}(v::SymplecticVector{n, T2}) where {n, T1<:Integer, T2<:Integer}
    return SymplecticVector{n, T1}(T1(v.a), T1(v.b))
end

convert(::Type{SymplecticVector{n, T}}, v::SymplecticVector{n, T}) where {n, T<:Integer} = v
convert(::Type{SymplecticVector{n, T1}}, v::SymplecticVector{n, T2}) where {n, T1<:Integer, T2<:Integer} = SymplecticVector{n, T1}(v)
promote_rule(::Type{SymplecticVector{n, T1}}, ::Type{SymplecticVector{n, T2}}) where {n, T1<:Integer, T2<:Integer} = SymplecticVector{n, promote_rule(T1, T2)}


############################################################################################
## Methods for extracting the data stored in a SymplecticVector{n, T} in different forms  ##
############################################################################################

eltype(::Type{SymplecticVector{n, T}}) where {n, T} = T
eltype(::SymplecticVector{n, T}) where {n, T} = T

"""
    halfdimension(v::SymplecticVector{n, T}) where {n, T} = n

For a vector ``v∈ℤ₂²ⁿ`` represented by a ``SymplecticVector{n, T}``, return ``n``.
"""
halfdimension(v::SymplecticVector{n, T}) where {n, T} = n

"""
    dimension(v::SymplecticVector{n, T}) where {n, T} = 2n

For a vector ``v∈ℤ₂²ⁿ`` represented by a SymplecticVector{n, T}``, return ``2n``.
"""
dimension(v::SymplecticVector{n, T}) where {n, T} = 2n

"""
    data(v::SymplecticVector{n, T})

Returns the internal representation of the vector ``v∈ℤ₂²ⁿ`` stored by `v`.

Objects of type `SymplecticVector{n, T}` represent a symplectic vector ``v∈ℤ₂²ⁿ`` through
the bits in the binary expansions of a pair of integers of type `T`.  These integers are
returned by data(v).

See also [`SymplecticVector{n, T}`](@ref).

# Examples
```jldoctest
julia> v = SymplecticVector{3, UInt8}(7, 4)
101011
julia> data(v)
(0x07, 0x04)
```
"""
function data(v::SymplecticVector{n, T}) where {n, T}
    return (v.a, v.b)
end

"""
    vector(v::SymplecticVector{n, T}) where {n, T}

Returns the symplectic vector ``v∈ℤ₂²ⁿ`` as a ``NTuple{2n, UInt8}``.

# Examples
```jldoctest
julia> v = SymplecticVector{3, UInt8}(7, 4)
101011
julia> vector(v)
(0x01, 0x00, 0x01, 0x00, 0x01, 0x01)
"""
function vector(v::SymplecticVector{n, T}) where {n, T}
    vec = []
    for k = 0:(n-1)
        push!(vec, (v.a>>>k) & 1)
        push!(vec, (v.b>>>k) & 1)
    end
    return NTuple{2n, UInt8}(vec)
end

"""
    bitstring(v::SymplecticVector{n, T}) where {n, T}

Returns the symplectic vector ``v∈ℤ₂²ⁿ`` as a string of bits.

# Examples
```jldoctest
julia> v = SymplecticVector{3, UInt8}(7, 4)
101011
julia> bitstring(v)
"101011"
```
"""
function bitstring(v::SymplecticVector{n, T}) where {n, T}
    abits = reverse(bitstring(v.a)[end-n+1:end])
    bbits = reverse(bitstring(v.b)[end-n+1:end])
    return join(Base.Iterators.flatten(zip([abits...], [bbits...])))
end

# More efficient bitstring computations for specific data types T. Interleavebits methods
# are defined in utils.jl
function bitstring(v::SymplecticVector{n, UInt8}) where {n}
    return reverse(bitstring(interleavebits(v.a, v.b))[end-2n+1:end])
end

function bitstring(v::SymplecticVector{n, UInt16}) where {n}
    return reverse(bitstring(interleavebits(v.a, v.b))[end-2n+1:end])
end

function bitstring(v::SymplecticVector{n, UInt32}) where {n}
    return reverse(bitstring(interleavebits(v.a, v.b))[end-2n+1:end])
end

function bitstring(v::SymplecticVector{n, UInt64}) where {n}
    return reverse(bitstring(interleavebits(v.a, v.b))[end-2n+1:end])
end

function bitstring(v::SymplecticVector{n, UInt128}) where {n}
    a1 = UInt64(v.a & 0x0000000000000000ffffffffffffffff)
    a2 = UInt64(v.a >> 64)
    b1 = UInt64(v.b & 0x0000000000000000ffffffffffffffff)
    b2 = UInt64(v.b >> 64)

    bits = bitstring(interleavebits(a2, b2)) * bitstring(interleavebits(a1, b1))
    return reverse(bits[end-2n+1:end])
end

# Overloading `bitstring` function for printing SymplecticVector{n,T} when T is BigInt
function bitstring(x::BigInt)
    bits = ""
    while x != 0
        bits = string(x & 1) * bits
        x >>>= 1
    end
    return bits
end

function bitstring(v::SymplecticVector{n, BigInt}) where {n}
    abits = rpad(reverse(bitstring(v.a)), n, '0')
    bbits = rpad(reverse(bitstring(v.b)), n, '0')
    return join(Base.Iterators.flatten(zip([abits...], [bbits...])))
end


############################################################################################
## Printing.                                                                              ##
## Although internally vectors are stored as integers, it's more convenient to output the ##
## bits in the binary expansion (i.e. the vector itself) rather than the integer.         ##
############################################################################################

function show(io::IO, v::SymplecticVector{n, T}) where {n, T}
    print(io, bitstring(v))
end


############################################################################################
## Functions for generating random symplectic vectors.                                    ##
############################################################################################

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, T}}) where {n, T}
    return SymplecticVector{n, T}(rand(rng, 0:big(2)^n-1, 2)...)
end

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, T}}, dims...) where {n, T}
    return SymplecticVector{n, T}.(
        rand(rng, 0:big(2)^n-1, dims...), rand(rng, 0:big(2)^n-1, dims...)
        )
end

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n}}) where {n}
    T = typerequired(n)
    return rand(SymplecticVector{n, T})
end

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n}}, dims...) where {n}
    T = typerequired(n)
    return rand(SymplecticVector{n, T}, dims...)
end


############################################################################################
## Functions defining basis arithmetic operations for vectors in ℤ₂²ⁿ.                    ##
############################################################################################

zero(::Type{SymplecticVector{n, T}}) where {n, T} = SymplecticVector{n, T}(zero(T), zero(T))

function zero(::Type{SymplecticVector{n}}) where {n}
    T = typerequired(n)
    return zero(SymplecticVector{n, T})
end

function iszero(u::SymplecticVector{n, T}) where {n, T}
    if iszero(u.a) && iszero(u.b)
        return true
    else
        return false
    end
end

function +(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return SymplecticVector{n, T}(u.a ⊻ v.a, u.b ⊻ v.b)
end

-(v::SymplecticVector{n, T}) where {n, T} = v
-(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T} = u + v

## In ℤ₂ᴺ, the only scalars are 0 and 1.
function *(k::Integer, v::SymplecticVector{n, T}) where {n, T}
    if k & 1 == 1
        return v
    end
    return zero(SymplecticVector{n, T})
end

*(v::SymplecticVector{n, T}, k::Integer) where {n, T} = *(k, v)

function dotproduct(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return parity((u.a & v.a) ⊻ (u.b & v.b))
end

⋅(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T} = dotproduct(u, v)

"""
    symplecticform(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}

Compute the symplectic inner product of the symplectic vectors ``u,v∈ℤ₂ⁿ×ℤ₂ⁿ``.

The symplectic form is defined as `sum((u[2k-1] && v[2k]) ⊻ (u[2k] && v[2k-1]) for k=1:n)`
(all arithmetic in ``ℤ₂``).

# Examples
```jldoctest
julia> u = SymplecticVector{3}(7, 4)
101011
julia> v = SymplecticVector{3}(5, 6)
100111
julia> symplecticform(u, v)
0x01
```
"""
function symplecticform(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return parity((u.a & v.b) ⊻ (u.b & v.a))
end

"""
    ⋆(u, v)

Alias for `symplecticform`, `u⋆v` calls `symplecticform(u, v)`.
"""
function ⋆(u, v)
    return symplecticform(u, v)
end

"""
    lift(u::SymplecticVector{n, T1}, m::Integer, T2) where {n, T1<:Integer}

Given an element u∈ℤ₂²ⁿ, return u⊕0∈ℤ₂²⁽ⁿ⁺ᵐ⁾.

The input is a `SymplecticVector{n, T1}` for some type `T1`, the output has type
`SymplectiCVector{n+m, T2}`.
"""
function lift(u::SymplecticVector{n, T1}, m::Integer, T2) where {n, T1<:Integer}
    return SymplecticVector{n+m, T2}(data(u)...)
end

function lift(u::SymplecticVector{n, T}, m::Integer) where {n, T<:Integer}
    T2 = typerequired(n+m)
    return SymplecticVector{n+m, T2}(data(u)...)
end

############################################################################################
## Comparison operators                                                                   ##
############################################################################################

function hash(u::SymplecticVector{n, T}) where {n, T<:Integer}
    return hash(data(u))
end

"""
    isequal(u, v)

Check if two symplectic vectors are equal (regardless of type parameters). E.g. 01 == 0100.

# Examples
```jldoctest
julia> u = SymplecticVector{3, UInt8}(7, 4)
101011
julia> v = SymplecticVector{3, UInt8}(5, 6)
100111
julia> isequal(u, v)
false
julia> v = SymplecticVector{9, UInt16}(7, 4)
101011000000000000
julia> isequal(u, v)
true
"""
function isequal(u::SymplecticVector{n1, T1}, v::SymplecticVector{n2, T2}) where {n1, n2, T1<:Integer, T2<:Integer}
    if isequal(u.a, v.a) && isequal(u.b, v.b)
        return true
    else
        return false
    end
end

==(u::SymplecticVector{n1, T1}, v::SymplecticVector{n2, T2}) where {n1, n2, T1<:Integer, T2<:Integer} = isequal(u, v)

"""
    isless(u, v)

Defines a partial order on the set of Symplectic vectors.

00,10,01,11,0010,1010,0110,1110,0001,1001,0101,1101,0011,1011,0111,1111,000010,100010
"""
function isless(u::SymplecticVector{n1, T1}, v::SymplecticVector{n2, T2}) where {n1, n2, T1<:Integer, T2<:Integer}
    if interleavebits(u.a, u.b) < interleavebits(v.a, v.b)
        return true
    else
        return false
    end
end


############################################################################################
## Symplectic group stuff                                                                 ##
############################################################################################

"""
    symplecticgrouporder(n)

Compute the order of the symplectic group ``Sp(2n,ℤ₂)``.
"""
function symplecticgrouporder(n)
    order = BigInt(2)^(n^2)
    for k in 1:n
        order *= BigInt(4)^k - 1
    end
    if order ≤ typemax(Int64)
        return Int64(order)
    elseif order ≤ typemax(Int128)
        return Int128(order)
    end
    return order
end


############################################################################################
## The symplectic Gram-Schmidt procedure                                                  ##
############################################################################################

"""
    symplecticgramschmidt(Ω::Vector{SymplecticVector{n, T}})

For as set of vectors ``Ω``, return a symplectic basis for ``span(Ω \\ (Ω ∩ Ω^⟂))``.
"""
function symplecticgramschmidt(Ω::Vector{T}) where {T<:AbstractSymplecticVector}
    basis::Vector{T} = []
    while length(Ω) ≠ 0
        symplecticpairflag = false
        u = Ω[1]
        for k = 2:length(Ω)
            v = Ω[k]
            if u ⋆ v == 1
                symplecticpairflag = true
                push!(basis, u)
                push!(basis, v)
                Ω = vcat(Ω[2:k-1], Ω[k+1:end])
                for k = 1:length(Ω)
                    w = Ω[k]
                    Ω[k] = w + ((v ⋆ w) * u) + ((u ⋆ w) * v)
                end
                break
            end
        end
        if !(symplecticpairflag)
            Ω = Ω[2:end]
        end
    end
    return basis
end


############################################################################################
## Types for symplectic maps                                                              ##
############################################################################################

abstract type AbstractSymplecticMap end


############################################################################################
## Symplectic Transvections                                                               ##
############################################################################################

"""
    Transvection{n, T}

A type that stores a symplectic transvection.

This is a map of the form ``Zₕ(v)=v+(h⋆v)h`` where ``h⋆v`` is the symplectic inner product
between ``h`` and ``v``. Internally, a transvection ``Zₕ`` of type `Transvection{n, T}`
stores ``h`` as a symplectic vector of type `SymplecticVector{n, T}`.

See also ['SymplecticVector{n, T}'](@ref).
"""
struct Transvection{n, T} <: AbstractSymplecticMap
    h::SymplecticVector{n, T}
    function Transvection{n, T}(h::SymplecticVector{n, T}) where {n, T<:Integer}
        return new(h)
    end
end

function Transvection{n}(h::SymplecticVector{n, T}) where {n, T<:Integer}
    return Transvection{n, T}(h)
end

function Transvection(h::SymplecticVector{n, T}) where {n, T<:Integer}
    return Transvection{n, T}(h)
end

## Symplectic transvections get printed like their underlying symplectic vectors
function show(io::IO, h::Transvection{n, T}) where {n, T}
    print(io, bitstring(h.h))
end

"""
    transvection(h, v)

Compute the symplectic transvection ``Zₕv:=v+[h,v]h.

If `h` is a vector then the transvections in `h` get applied sequentially. If
`v` is a vector then the transvection(s) gets applied to `v` element-wise.

See also [`findtransvection`](@ref).

# Examples
```jldoctest
julia> a = SymplecticVector{3}(11)
110100
julia> b = SymplecticVector{3}(17)
010110
julia> h = Transvection{3}(a)
110100
julia> transvection(h, b)
100010
```
"""
function transvection(h, v) end

function transvection(h::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    if h ⋆ v == zero(T)
        return v
    else
        return h + v
    end
end

function transvection(h::Transvection{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return transvection(h.h, v)
end

function transvection(H::Vector, v::SymplecticVector{n, T}) where {n, T}
    for h in H
        v = transvection(h, v)
    end
    return v
end

function transvection(h::Transvection{n, T}, V::Vector) where {n, T}
    return [transvection(h, v) for v in V]
end

function transvection(H::Vector, V::Vector)
    return [transvection(H, v) for v in V]
end

*(h::Transvection{n, T}, v::SymplecticVector{n, T}) where {n, T} = transvection(h, v)


## Random transvections
function rand(rng::AbstractRNG, ::SamplerType{Transvection{n, T}}) where {n, T}
    return Transvection{n, T}(rand(SymplecticVector{n, T}))
end

function rand(rng::AbstractRNG, ::SamplerType{Transvection{n, T}}, dims...) where {n, T}
    return Transvection{n, T}.(rand(SymplecticVector{n, T}, dims...))
end

function rand(rng::AbstractRNG, ::SamplerType{Transvection{n}}) where {n}
    T = typerequired(n)
    return rand(Transvection{n, T})
end

function rand(rng::AbstractRNG, ::SamplerType{Transvection{n}}, dims...) where {n}
    T = typerequired(n)
    return rand(Transvection{n, T}, dims...)
end

"""
    findtransvection(u::SymplecticVector{n, T}, v::SymplecticVector{n, T})

Find vectors ``s,t∈ℤ₂²ⁿ`` such that ``v=ZₛZₜu``.

Returns a vector `[s, t]` where `s` and `t` are of type `Transvection{n, T}`. This procedure
for finding transvections with the required properties is described in the proof of Lemma 2
in J. Math. Phys. 55, 122202 (2014). If ``u⋆v=1`` then ``t=0``, if ``u=v`` then ``s=0=t``.

See also [`Transvection{n, T}`](@ref).
"""
function findtransvection(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    if u == v
        return repeat([Transvection{n, T}(SymplecticVector{n, T}(zero(T), zero(T)))], 2)
    elseif symplecticform(u, v) == one(T)
        return [Transvection{n, T}(u+v), Transvection{n, T}(SymplecticVector{n, T}(zero(T), zero(T)))]
    end
    for j=0:(n-1)
        u₁ = (u.a >>> j) & 1; u₂ = (u.b >>> j) & 1
        v₁ = (v.a >>> j) & 1; v₂ = (v.b >>> j) & 1
        if ((u₁ | u₂) & (v₁ | v₂)) == 1
            for q = 0:3
                if (u₁&q ⊻ u₂&(q>>>1)) == 1 == (v₁&q ⊻ v₂&(q>>>1))
                    w::SymplecticVector{n, T} = u + SymplecticVector{n, T}((q>>>1)<<j, (q&1)<<j)
                    return [Transvection{n, T}(u+w), Transvection{n, T}(v+w)]
                end
            end
        end
    end
    for j=0:(n-1)
        u₁ = (u.a >>> j) & 1; u₂ = (u.b >>> j) & 1
        if (u₁ | u₂) == 1
            for k = 0:(n-1)
                v₁ = (v.a >>> k) & 1; v₂ = (v.b >>> k) & 1
                if (v₁ | v₂) == 1
                    for q = 0:15
                        if (u₁&q ⊻ u₂&(q>>>1)) == 1 == (v₁&(q>>>2) ⊻ v₂&(q>>>3))
                            w::SymplecticVector{n, T} = u + SymplecticVector{n, T}((((q>>>1)&1)<<j) | (((q>>>3)&1)<<k), ((q&1)<<j) | (((q>>>2)&1)<<k))
                            return [Transvection{n, T}(u+w), Transvection{n, T}(v+w)]
                        end
                    end
                end
            end
        end
    end
end


############################################################################################
## Subspaces                                                                              ##
############################################################################################

struct Subspace{n, dim, T}
    basis::NTuple{dim, SymplecticVector{n, T}}
    function Subspace{n, dim, T}(b::NTuple{dim, SymplecticVector{n, T}}) where {n, dim, T<:Integer}
        return new(b)
    end
end

function dimension(S::Subspace{n, dim, T}) where {n, dim, T<:Integer}
    return dim
end

function isisotropic(S::Subspace{n, dim, T}) where {n, dim, T<:Integer}
    for i = 1:dim
        for j = 1:dim
            if !iszero(S.basis[i] ⋆ S.basis[j])
                return false
            end
        end
    end
    return true
end

function isLangrangian(S::Subspace{n, dim, T}) where {n, dim, T<:Integer}
    if dim == n && isisotropic(S)
        return true
    else
        return false
    end
end

"""
    SYMPLECTICImproved(n, i)

The SYMPLECTICImproved algorithm from J. Math. Phys. 55, 122202 (2014).
"""
function SYMPLECTICImproved(n, i)
    type = typerequired(n)
    s = big(1) << 2n - 1
    k = i % s + 1

    f₁ = SymplecticVector{n, type}(deinterleavebits(k)...)
    
    T = findtransvection(SymplecticVector{n, type}(one(type),zero(type)), f₁)

    b = swapbits(BigInt(floor(i / s)) << 1, 0, 1) & (big(1) << 2n - 1)

    e = SymplecticVector{n, type}(reverse(deinterleavebits(b))...)
    h₀ = transvection(T, e)

    if b & 1 == 1
        T¹ = [Transvection{n}(h₀)]
    else
        T¹ = [Transvection{n}(h₀), Transvection{n}(f₁)]
    end
    f₂ = transvection([T, T¹...], SymplecticVector{n, type}(zero(type),one(type)))

    if n == 1
        return [f₁, f₂]
    else
        return transvection([T, T¹...],
                [
                    [lift(u, 1) for u in SYMPLECTICImproved(n-1, BigInt(floor(i / s)) >> (2n-1))]...,
                    SymplecticVector{n, type}(one(type) << (n-1), 0),
                    SymplecticVector{n, type}(0, one(type) << (n-1))
                ]
            )
    end
end

end  # module
