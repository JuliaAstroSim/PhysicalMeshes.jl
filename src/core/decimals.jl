#! This file is deprecated

# Quantity
decimal(x::Quantity{T,D,U}) where {T<:FloatTypes, D, U} = 
    Quantity{Decimal, D, U}(x.val)

floatnumber(x::Number) = x

# Always return Float64
function floatnumber(x::Decimal)
    str = string(x)
    return parse(Float64, str)
end

floatnumber(x::Quantity{T,D,U}) where {T<:Decimal, D, U} = 
    Quantity{Float64, D, U}(floatnumber(x.val))

# PVector
decimal(p::PVector2D) = PVector2D(decimal(p.x), decimal(p.y))

decimal(p::PVector) = PVector(decimal(p.x), decimal(p.y), decimal(p.z))

floatnumber(p::PVector2D) = PVector2D(floatnumber(p.x), floatnumber(p.y))

floatnumber(p::PVector) = PVector(floatnumber(p.x), floatnumber(p.y), floatnumber(p.z))