# Mostly bit twiddling hacks from
# http://www-graphics.stanford.edu/~seander/bithacks.html?fbclid=IwAR1xlKWawUWqNKpAZgWNN9-1FLPCXCnhx6f88h7_OUH4Js0wzxztwcEqpJY

function parity(a::Integer)
    parity = 0
    while a ≠ 0
        parity ⊻= 1
        a = a & (a - 1)
    end
    return parity
end

function parity(a::UInt8)
    a = (a ⊻ (a >> 4)) & 0xf
    return (0x6996 >> a) & 0x1
end

function parity(a::UInt16)
    a ⊻= a >> 8
    a ⊻= a >> 4
    a &= 0xf
    return (0x6996 >> a) & 0x1
end

function parity(a::UInt32)
    a ⊻= a >> 16
    a ⊻= a >> 8
    a ⊻= a >> 4
    a &= 0xf
    return (0x6996 >> a) & 0x1
end

function parity(a::UInt64)
    a ⊻= a >> 32
    a ⊻= a >> 16
    a ⊻= a >> 8
    a ⊻= a >> 4
    a &= 0xf
    return (0x6996 >> a) & 0x1
end

function parity(a::UInt128)
    a ⊻= a >> 64
    a ⊻= a >> 32
    a ⊻= a >> 16
    a ⊻= a >> 8
    a ⊻= a >> 4
    a &= 0xf
    return (0x6996 >> a) & 0x1
end

function parity(a::BigInt)
    k = 0
    
    while a ≠ 0
        k ⊻= a & 1
        k >> 1
    end

    return k
end

function interleavebits(a::UInt8, b::UInt8)
    a = UInt16(a); b = UInt16(b);

    a = (a | (a << 4)) & 0x0f0f
    a = (a | (a << 2)) & 0x3333
    a = (a | (a << 1)) & 0x5555

    b = (b | (b << 4)) & 0x0f0f
    b = (b | (b << 2)) & 0x3333
    b = (b | (b << 1)) & 0x5555

    return a | (b << 1)
end

function interleavebits(a::UInt16, b::UInt16)
    a = UInt32(a); b = UInt32(b);

    a = (a | (a << 8)) & 0x00ff00ff
    a = (a | (a << 4)) & 0x0f0f0f0f
    a = (a | (a << 2)) & 0x33333333
    a = (a | (a << 1)) & 0x55555555

    b = (b | (b << 8)) & 0x00ff00ff
    b = (b | (b << 4)) & 0x0f0f0f0f
    b = (b | (b << 2)) & 0x33333333
    b = (b | (b << 1)) & 0x55555555

    return a | (b << 1)
end

function interleavebits(a::UInt32, b::UInt32)
    a = UInt64(a); b = UInt64(b);

    a = (a | (a << 16)) & 0x0000ffff0000ffff
    a = (a | (a <<  8)) & 0x00ff00ff00ff00ff
    a = (a | (a <<  4)) & 0x0f0f0f0f0f0f0f0f
    a = (a | (a <<  2)) & 0x3333333333333333
    a = (a | (a <<  1)) & 0x5555555555555555

    b = (b | (b << 16)) & 0x0000ffff0000ffff
    b = (b | (b <<  8)) & 0x00ff00ff00ff00ff
    b = (b | (b <<  4)) & 0x0f0f0f0f0f0f0f0f
    b = (b | (b <<  2)) & 0x3333333333333333
    b = (b | (b <<  1)) & 0x5555555555555555

    return a | (b << 1)
end

function interleavebits(a::UInt64, b::UInt64)
    a = UInt128(a); b = UInt128(b)

    a = (a | (a << 32)) & 0x00000000ffffffff00000000ffffffff
    a = (a | (a << 16)) & 0x0000ffff0000ffff0000ffff0000ffff
    a = (a | (a <<  8)) & 0x00ff00ff00ff00ff00ff00ff00ff00ff
    a = (a | (a <<  4)) & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
    a = (a | (a <<  2)) & 0x33333333333333333333333333333333
    a = (a | (a <<  1)) & 0x55555555555555555555555555555555

    b = (b | (b << 32)) & 0x00000000ffffffff00000000ffffffff
    b = (b | (b << 16)) & 0x0000ffff0000ffff0000ffff0000ffff
    b = (b | (b <<  8)) & 0x00ff00ff00ff00ff00ff00ff00ff00ff
    b = (b | (b <<  4)) & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
    b = (b | (b <<  2)) & 0x33333333333333333333333333333333
    b = (b | (b <<  1)) & 0x55555555555555555555555555555555

    return a | (b << 1)
end

function interleavebits(a::UInt128, b::UInt128)
    a = BigInt(a); b = BigInt(b)

    a = (a | (a << 64)) & ((BigInt(0x0000000000000000ffffffffffffffff) << 128) | (0x0000000000000000ffffffffffffffff))
    a = (a | (a << 32)) & ((BigInt(0x00000000ffffffff00000000ffffffff) << 128) | (0x00000000ffffffff00000000ffffffff))
    a = (a | (a << 16)) & ((BigInt(0x0000ffff0000ffff0000ffff0000ffff) << 128) | (0x0000ffff0000ffff0000ffff0000ffff))
    a = (a | (a <<  8)) & ((BigInt(0x00ff00ff00ff00ff00ff00ff00ff00ff) << 128) | (0x00ff00ff00ff00ff00ff00ff00ff00ff))
    a = (a | (a <<  4)) & ((BigInt(0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f) << 128) | (0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f))
    a = (a | (a <<  2)) & ((BigInt(0x33333333333333333333333333333333) << 128) | (0x33333333333333333333333333333333))
    a = (a | (a <<  1)) & ((BigInt(0x55555555555555555555555555555555) << 128) | (0x55555555555555555555555555555555))

    b = (b | (b << 64)) & ((BigInt(0x0000000000000000ffffffffffffffff) << 128) | (0x0000000000000000ffffffffffffffff))
    b = (b | (b << 32)) & ((BigInt(0x00000000ffffffff00000000ffffffff) << 128) | (0x00000000ffffffff00000000ffffffff))
    b = (b | (b << 16)) & ((BigInt(0x0000ffff0000ffff0000ffff0000ffff) << 128) | (0x0000ffff0000ffff0000ffff0000ffff))
    b = (b | (b <<  8)) & ((BigInt(0x00ff00ff00ff00ff00ff00ff00ff00ff) << 128) | (0x00ff00ff00ff00ff00ff00ff00ff00ff))
    b = (b | (b <<  4)) & ((BigInt(0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f) << 128) | (0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f))
    b = (b | (b <<  2)) & ((BigInt(0x33333333333333333333333333333333) << 128) | (0x33333333333333333333333333333333))
    b = (b | (b <<  1)) & ((BigInt(0x55555555555555555555555555555555) << 128) | (0x55555555555555555555555555555555))

    return a | (b << 1)
end

function interleavebits(a::BigInt, b::BigInt)
    c = zero(BigInt)
    k = 0

    while (a ≠ 0) || (b ≠ 0)
        c |= (a & 1) << (2k)
        c |= (b & 1) << (2k+1)

        a >>= 1
        b >>= 1
        k += 1
    end
    
    return c
end

function deinterleavebits(c::UInt8)
    a = c & 0x55; b = (c & 0xaa) >> 1;

    a = (a | (a >> 1)) & 0x33
    a = (a | (a >> 2)) & 0x0f

    b = (b | (b >> 1)) & 0x33
    b = (b | (b >> 2)) & 0x0f

    return (a, b)
end

function deinterleavebits(c::UInt16)
    a = c & 0x5555; b = (c & 0xaaaa) >> 1;

    a = (a | (a >> 1)) & 0x3333
    a = (a | (a >> 2)) & 0x0f0f
    a = (a | (a >> 4)) & 0x00ff

    b = (b | (b >> 1)) & 0x3333
    b = (b | (b >> 2)) & 0x0f0f
    b = (b | (b >> 4)) & 0x00ff

    return (UInt8(a), UInt8(b))
end

function deinterleavebits(c::UInt32)
    a = c & 0x55555555; b = (c & 0xaaaaaaaa) >> 1;

    a = (a | (a >> 1)) & 0x33333333
    a = (a | (a >> 2)) & 0x0f0f0f0f
    a = (a | (a >> 4)) & 0x00ff00ff
    a = (a | (a >> 8)) & 0x0000ffff

    b = (b | (b >> 1)) & 0x33333333
    b = (b | (b >> 2)) & 0x0f0f0f0f
    b = (b | (b >> 4)) & 0x00ff00ff
    b = (b | (b >> 8)) & 0x0000ffff

    return (UInt16(a), UInt16(b))
end

function deinterleavebits(c::UInt64)
    a = c & 0x5555555555555555; b = (c & 0xaaaaaaaaaaaaaaaa) >> 1;

    a = (a | (a >>  1)) & 0x3333333333333333
    a = (a | (a >>  2)) & 0x0f0f0f0f0f0f0f0f
    a = (a | (a >>  4)) & 0x00ff00ff00ff00ff
    a = (a | (a >>  8)) & 0x0000ffff0000ffff
    a = (a | (a >> 16)) & 0x00000000ffffffff

    b = (b | (b >>  1)) & 0x3333333333333333
    b = (b | (b >>  2)) & 0x0f0f0f0f0f0f0f0f
    b = (b | (b >>  4)) & 0x00ff00ff00ff00ff
    b = (b | (b >>  8)) & 0x0000ffff0000ffff
    b = (b | (b >> 16)) & 0x00000000ffffffff

    return (UInt32(a), UInt32(b))
end

function deinterleavebits(c::UInt128)
    a = c & 0x55555555555555555555555555555555; b = (c & 0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa) >> 1;

    a = (a | (a >>  1)) & 0x33333333333333333333333333333333
    a = (a | (a >>  2)) & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
    a = (a | (a >>  4)) & 0x00ff00ff00ff00ff00ff00ff00ff00ff
    a = (a | (a >>  8)) & 0x0000ffff0000ffff0000ffff0000ffff
    a = (a | (a >> 16)) & 0x00000000ffffffff00000000ffffffff
    a = (a | (a >> 32)) & 0x0000000000000000ffffffffffffffff

    b = (b | (b >>  1)) & 0x33333333333333333333333333333333
    b = (b | (b >>  2)) & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f
    b = (b | (b >>  4)) & 0x00ff00ff00ff00ff00ff00ff00ff00ff
    b = (b | (b >>  8)) & 0x0000ffff0000ffff0000ffff0000ffff
    b = (b | (b >> 16)) & 0x00000000ffffffff00000000ffffffff
    b = (b | (b >> 32)) & 0x0000000000000000ffffffffffffffff

    return (UInt64(a), UInt64(b))
end

function deinterleavebits(c::BigInt)
    a = zero(BigInt); b = zero(BigInt)
    k = 0

    while c ≠ 0
        a |= (c & 1) << k
        c >>= 1
        b |= (c & 1) << k
        c >>= 1

        k += 1
    end
    
    return (a, b)
end

function reversebits(a::UInt8)
    a = ((a >> 1) & 0x55) | ((a & 0x55) << 1)
    a = ((a >> 2) & 0x33) | ((a & 0x33) << 2)
    a = ((a >> 4)       ) | ((a       ) << 4)
    return a
end

function reversebits(a::UInt16)
    a = ((a >> 1) & 0x5555) | ((a & 0x5555) << 1)
    a = ((a >> 2) & 0x3333) | ((a & 0x3333) << 2)
    a = ((a >> 4) & 0x0f0f) | ((a & 0x0f0f) << 4)
    a = ((a >> 8)         ) | ((a         ) << 8)
    return a
end

function reversebits(a::UInt32)
    a = ((a >> 1 ) & 0x55555555) | ((a & 0x55555555) << 1 )
    a = ((a >> 2 ) & 0x33333333) | ((a & 0x33333333) << 2 )
    a = ((a >> 4 ) & 0x0f0f0f0f) | ((a & 0x0f0f0f0f) << 4 )
    a = ((a >> 8 ) & 0x00ff00ff) | ((a & 0x00ff00ff) << 8 )
    a = ((a >> 16)            ) | ((a             ) << 16)
    return a
end

function reversebits(a::UInt64)
    a = ((a >>  1) & 0x5555555555555555) | ((a & 0x5555555555555555) <<  1)
    a = ((a >>  2) & 0x3333333333333333) | ((a & 0x3333333333333333) <<  2)
    a = ((a >>  4) & 0x0f0f0f0f0f0f0f0f) | ((a & 0x0f0f0f0f0f0f0f0f) <<  4)
    a = ((a >>  8) & 0x00ff00ff00ff00ff) | ((a & 0x00ff00ff00ff00ff) <<  8)
    a = ((a >> 16) & 0x0000ffff0000ffff) | ((a & 0x0000ffff0000ffff) << 16)
    a = ((a >> 32)                     ) | ((a                     ) << 32)
    return a
end

function reversebits(a::UInt128)
    a = ((a >>  1) & 0x55555555555555555555555555555555) | ((a & 0x55555555555555555555555555555555) <<  1)
    a = ((a >>  2) & 0x33333333333333333333333333333333) | ((a & 0x33333333333333333333333333333333) <<  2)
    a = ((a >>  4) & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f) | ((a & 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f) <<  4)
    a = ((a >>  8) & 0x00ff00ff00ff00ff00ff00ff00ff00ff) | ((a & 0x00ff00ff00ff00ff00ff00ff00ff00ff) <<  8)
    a = ((a >> 16) & 0x0000ffff0000ffff0000ffff0000ffff) | ((a & 0x0000ffff0000ffff0000ffff0000ffff) << 16)
    a = ((a >> 32) & 0x00000000ffffffff00000000ffffffff) | ((a & 0x00000000ffffffff00000000ffffffff) << 32)
    a = ((a >> 64)                                     ) | ((a                                     ) << 64)
    return a
end

function reversebits(a::BigInt)
    maxbit = Integer(floor(log2(a)))
    b = zero(BigInt)

    for n=0:maxbit
        b |= ((a >> n) & 1) << (maxbit - n)
    end

    return b
end

function swapbits(a, i, j)
    x = ((a >> i) ⊻ (a >> j)) & ((one(a) << 1) - one(a))
    r = a ⊻ ((x << i) | (x << j));
    return r
end