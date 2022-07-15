"""
    BinarySymplectic.jl

Tools for working with symplectic vector spaces and symplectic groups over ℤ₂.

https://github.com/mzurel/BinarySymplectic.jl
"""
module BinarySymplectic

import Base: (==), (+), (-), (*), (/), (^)
import Base: convert, promote, eltype, bitstring, show, hash, rand, zero
import Random: AbstractRNG, SamplerType

export SymplecticVector, SymplecticMap
export halfdimension, dimension, data, vector, bitstring
export symplecticform, dotproduct
export symplecticmap, symplectictransvection
export symplecticgrouporder

include("utils.jl")

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
## TODO: can optimize here by creating a distinct algorithm for array generation.         ##
############################################################################################
function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, T}}) where {n, T}
    return SymplecticVector{n, T}(rand(rng, 0:big(2)^n-1, 2)...)
end

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n}}) where {n}
    return rand(SymplecticVector{n, typerequired(n)})
end


############################################################################################
## Functions defining basis arithmetic operations for vectors in ℤ₂²ⁿ.                    ##
############################################################################################
zero(::Type{SymplecticVector{n, T}}) where {n, T} = SymplecticVector{n, T}(zero(T), zero(T))
function zero(::Type{SymplecticVector{n}}) where {n}
    T = typerequired(n)
    return zero(SymplecticVector{n, T})
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


############################
## Symplectic group stuff ##
############################
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

end  # module
