# Quantity
decimal(x::Quantity{T,D,U}) where {T<:FloatTypes, D, U} = 
    Quantity{Decimal, D, U}(x.val)

number(x::Number) = x

# Always return Float64
function number(x::Decimal)
    str = string(x)
    return parse(Float64, str)
end

number(x::Quantity{T,D,U}) where {T<:Decimal, D, U} = 
    Quantity{Float64, D, U}(number(x.val))

# PVector
decimal(p::PVector2D) = PVector2D(decimal(p.x), decimal(p.y))

decimal(p::PVector) = PVector(decimal(p.x), decimal(p.y), decimal(p.z))

number(p::PVector2D) = PVector2D(number(p.x), number(p.y))

number(p::PVector) = PVector(number(p.x), number(p.y), number(p.z))