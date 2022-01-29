"""
    BinarySymplectic

Tools for working with symplectic vector spaces and symplectic groups over ℤ₂.

https://github.com/mzurel/BinarySymplectic.jl
"""
module BinarySymplectic

import Base: (==), (+), (-), (*), (inv), (/), (//), (^), hash, show, rand, bitstring

export typerequired, SymplecticVector, vector, symplecticform, transvection, findtransvection, symplecticgrouporder, symplectic

"""
    typerequired(n::Integer)

Compute the smallest UInt type required to represent an element of ``ℤ₂²ⁿ``.

Elements of the symplectic vector space ``ℤ₂²ⁿ`` are represented as the bits in the binary
expansion of an integer. This function computes the smallest unsigned integer type with
enough bits to store an element of ``ℤ₂²ⁿ``. If ``n>64``, built-in unsigned types are too
small so BigInt is returned.

TODO: see if larger fixed width unsigned types from BitIntegers.jl will work for ``n>64``.
"""
function typerequired(n::Integer)
    if n ≤ 64
        bitsrequired = max(8, convert(Integer, 2^ceil(log2(2 * n))))
        return eval(Meta.parse("UInt"*string(bitsrequired)))
    else
        return BigInt
    end
end

"""
    SymplecticVector{n, T}

A type that stores an element of the symplectic vector space ``ℤ₂²ⁿ``.

The vector is stored as the bits in the binary expansion of an integer of type `T`. An
integer `a` is mapped to ``b₁e₁+b₂f₁+b₃e₂+b₄f₂+...`` where ``e₁``,``f₁``,``e₂``,``f₂``,...
are the standard symplectic basis vectors and ``b₁`` is the least significant bit in the
binary expansion of `a`, ``b₂`` is the second least significant bit, and so on.

# Examples
```jldoctest
julia> a = SymplecticVector{3, UInt8}(11)
110100
```
"""
struct SymplecticVector{n, T}
    vector::T
    function SymplecticVector{n, T}(vector::T) where {n, T<:Integer}
        return new(vector)
    end
end

function SymplecticVector{n, T1}(vector::T2) where {n, T1<:Integer, T2<:Integer}
    return SymplecticVector{n, T1}(convert(T1, vector))
end

"""
    SymplecticVector{n}(vector::Integer)

If the type `T` is not specified then the smallest unsigned type with enough bits to store
the vector is used.

See also [`typerequired`](@ref).

# Examples
```jldoctest
julia> a = SymplecticVector{3}(11)
110100
```
"""
@inline function SymplecticVector{n}(vector::Integer) where {n}
    T = typerequired(n)
    return SymplecticVector{n, T}(vector)
end

"""
    vector(v::SymplecticVector)

Returns the field `vector` of a object of type `SymplecticVector{n, T}`.

The field `vector` is an integer of type `T` representing the symplectic vector ``v∈ℤ₂²ⁿ``
through its binary expansion.

See also [`SymplecticVector{n, T}`](@ref).

# Examples
```jldoctest
julia> a = SymplecticVector{3, UInt8}(11)
110100
julia> vector(a)
0x0b
```
"""
function vector(v::SymplecticVector)
    return v.vector
end


############################################################################################
## Printing.                                                                              ##
## Although internally vectors are stored as integers, it's more convenient to output the ##
## bits in the binary expansion (i.e. the vector itself) rather than the integer.         ##
############################################################################################
function show(io::IO, v::SymplecticVector{n, T}) where {n, T}
    print(io, reverse(bitstring(v.vector)[end-(2*n)+1:end]))
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

function show(io::IO, v::SymplecticVector{n, BigInt}) where {n}
    print(io, reverse(lpad(bitstring(v.vector), 2*n, "0")))
end


############################################################################################
## Functions for generating random symplectic vectors.                                    ##
############################################################################################
rand(::Type{SymplecticVector{n, T}}) where {n, T} = SymplecticVector{n, T}(rand(0:big(2)^(2*n)-1))
rand(::Type{SymplecticVector{n}}) where {n} = rand(SymplecticVector{n, typerequired(n)})


############################################################################################
## Functions defining basis arithmetic operations for vectors in ℤ₂\^2n.                  ##
############################################################################################
function +(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return SymplecticVector{n, T}(u.vector ⊻ v.vector)
end

function *(k::Integer, v::SymplecticVector{n, T}) where {n, T}
    return SymplecticVector{n, T}(k .* v.vector)
end

*(v::SymplecticVector{n, T}, k::Integer) where {n, T} = *(k, v)


"""
    symplecticform(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}

Compute the symplectic inner product of the symplectic vectors ``u,v∈ℤ₂²ⁿ``.

The symplectic form is defined as `sum((u[2k-1] && v[2k]) ⊻ (u[2k] && v[2k-1]) for k=1:n)`
(all arithmetic in ``ℤ₂``).

# Examples
```jldoctest
julia> a = SymplecticVector{3}(11)
110100
julia> b = SymplecticVector{3}(17)
100010
julia> symplecticform(a, b)
0x01
```
"""
function symplecticform(u, v) end

"""
    symplecticform(n::Integer, u::T, v::T) where {T<:Integer}

If `u`, `v` are integers rather than `SymplecticVector{n, T}` types, then an extra argument
`n` is required to specify the dimension of the symplectic vector space ℤ₂²ⁿ.

# Examples
```jldoctest
julia> a = 11
11
julia> b = 17
17
julia> symplecticform(3, a, b)
1
```
"""
function symplecticform(n::Integer, u::T, v::T)::T where {T<:Integer}
    return reduce(⊻, ((u >>> (2*k-2)) & (v >>> (2*k-1)) & 1) ⊻ ((u >>> (2*k-1)) & (v >>> (2*k-2)) & 1) for k = 1:n)
end

function symplecticform(u::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return symplecticform(n, u.vector, v.vector)
end

"""
    ⋆(u, v)

Alias for `symplecticform`, `u⋆v` calls `symplecticform(u, v)`.
"""
function ⋆(u, v)
    return symplecticform(u, v)
end


###################################
## Symplectic transvection stuff ##
###################################
"""
    transvection(h, v)

Compute the symplectic transvection ``Zₕv:=v+[h,v]h``.

See also [`findtransvection`](@ref).

# Examples
```jldoctest
julia> a = SymplecticVector{3}(11)
110100
julia> b = SymplecticVector{3}(17)
010110
julia> transvection(a, b)
010110
```
"""
function transvection(h, v) end

function transvection(h::SymplecticVector{n, T}, v::SymplecticVector{n, T}) where {n, T}
    return v + ((h ⋆ v) * h)
end

"""
    transvection(H::Vector, v)

If the first argument is a `Vector` then the transvections labeled by elements of `H` are
applied to `v` sequentially.  I.e. the result is ``...Z_h₂ Z_h₁ v`` where `H=[h₁,h₂,...]`.
"""
function transvection(h::Vector{SymplecticVector{n, T}}, v::SymplecticVector{n, T}) where {n, T}
    for k in h
        v = transvection(k, v)
    end
    return v
end

"""
    transvection(h, V::Vector)

If the second argument is a vector then the transvection(s) `h` is(are) applied to each
element of `V`.
"""
transvection(h,V::Vector{SymplecticVector{n,T}}) where {n, T}=[transvection(h,v) for v in V]

"""
    transvection(n, h, v)

When `h` and `v` are integers rather than `SymplecticVector{n, T}` types, an extra argument
`n` is required to specify the vector space ℤ₂²ⁿ where the vectors `h` and `v` live. In this
case, `h` and `v` can also be vectors of integers with the same behaviour as before.
"""
function transvection(n, h, v) end

function transvection(n::Integer, h::T, v::T)::T where {T<:Integer}
    return v ⊻ (symplecticform(n, h, v) * h)
end

function transvection(n::Integer, H::Vector{T}, v::T)::T where {T<:Integer}
    for h in H
        v = transvection(n, h, v)
    end
    return v
end

function transvection(n::Integer, h, v::Vector{T})::Vector{T} where {T<:Integer}
    return [transvection(n, h, u) for u in v]
end


"""
    findtransvection(u::SymplecticVector{n, T}, v::SymplecticVector{n, T})

Find vectors h₁,h₂ ∈ ℤ₂²ⁿ such that v=Z_h₁ Z_h₂ u.

This procedure for finding transvections with the required properties is described in the
proof of Lemma 2 in J. Math. Phys. 55, 122202 (2014).

See also [`transvection`](@ref).
"""
function findtransvection(u, v) end


"""
    findtransvection(n::Integer, x::T, y::T) where {T<:Integer}

When `x` and `y` are integers rather than `SymplecticVector{n, T}` types, an extra argument
`n` is required to specify the vector space ℤ₂²ⁿ where the vectors `x` and `y` live.
"""
function findtransvection(n::Integer, x::T, y::T)::Vector{T} where {T<:Integer}
    if x == y
        return [zero(T), zero(T)]
    elseif symplecticform(n, x, y) == one(T)
        return [x⊻y, zero(T)]
    end
    for j in 1:n
        x1 = (x >>> (2*j-2)) & 1; x2 = (x >>> (2*j-1)) & 1
        y1 = (y >>> (2*j-2)) & 1; y2 = (y >>> (2*j-1)) & 1
        if ((x1 | x2) & (y1 | y2)) == 1
            for v in 0:3
                if (x1&v ⊻ x2&(v>>>1) == 1) && (y1&v ⊻ y2&(v>>>1) == 1)
                    z::T = x ⊻ ((v>>>1)<<(2*j-2)) ⊻ ((v&1)<<(2*j-1))
                    return [x⊻z, z⊻y]
                end
            end
        end
    end
    for j in 1:n
        x1 = (x >>> (2*j-2)) & 1; x2 = (x >>> (2*j-1)) & 1
        if (x1 | x2) == 1
            for k in 1:n
                y1 = (y >>> (2*k-2)) & 1; y2 = (y >>> (2*k-1)) & 1
                if (y1 | y2) == 1
                    for v in 0:15
                        if (x1&v ⊻ x2&(v>>>1) == 1) && (y1&(v>>>2) ⊻ y2&(v>>>3) == 1)
                            z::T = x ⊻ (((v>>>1)&1)<<(2*j-2)) ⊻ ((v&1)<<(2*j-1)) ⊻ ((v>>>3)<<(2*k-2)) ⊻ (((v>>>2)&1)<<(2*k-1))
                            return [x⊻z, z⊻y]
                        end
                    end
                end
            end
        end
    end
end

function findtransvection(u::SymplecticVector{n, T}, v::SymplecticVector{n, T})::Vector{SymplecticVector{n, T}} where {n, T}
    SymplecticVector{n, T}.(findtransvection(n, u.vector, v.vector))
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

"""
    SYMPLECTICImproved(nQubits, i)

Return an element of the symplectic group Sp(2n,ℤ₂) uniquely identified by ``i∈1:|Sp(2n,ℤ₂)|``.

The output is a vector of integers `V` representing binary symplectic vectors through their
binary expansions.  The symplectic group element is defined by the map ``e₁->V[1]``,
``f₁->V[2]``, ``e₂->V[3]``,...  The unique map from integers ``i∈1:|Sp(2n,ℤ₂)|`` is derived
from the SYMPLECTICImproved algorithm from J. Math. Phys. 55, 122202 (2014).
"""
function SYMPLECTICImproved(n::Integer, i::Integer)
    T = typerequired(n)
    s::T = big(2)^(2 * n) - 1
    k::T = (i % s) + 1

    e1 = one(T)
    t = findtransvection(n, e1, k)

    b::T = floor(T, i / s)

    e::T = e1 ⊻ (b>>>1)<<2
    h0 = transvection(n, t, e)

    if (b & 1) == 1
        tprime = [h0...]
    else
        tprime = [h0...,k]
    end

    if n == 1
        return [k, transvection(n, [t; tprime], T(2))]  # [f1,f2]
    else
        return transvection(n, [t; tprime], [one(T)<<(2*n-1), one(T)<<(2*n-2), SYMPLECTICImproved(n-1, T(floor(i / s)) >>> (2*n - 1))...])
    end
end

"""
    symplectic(n::Integer, i::Integer)

Return an element of the symplectic group Sp(2n,ℤ₂) uniquely identified by ``i∈1:|Sp(2n,ℤ₂)|``.

The output is a vector `V` of objects of type `SymplecticVector{n, T}`. The symplectic group
element is defined by the map ``e₁->V[1]``, ``f₁->V[2]``, ``e₂->V[3]``,....  The unique map
from integers ``i∈1:|Sp(2n,ℤ₂)|`` is derived from the SYMPLECTICImproved algorithm from
J. Math. Phys. 55, 122202 (2014).
"""
function symplectic(n::Integer, i::Integer)
    T = typerequired(n)
    return SymplecticVector{n, T}.(SYMPLECTICImproved(n, i))
end

end  # module
