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

    a | (b << 1)
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
