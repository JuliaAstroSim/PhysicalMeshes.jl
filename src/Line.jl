struct Line2D{T<:Number} <: AbstractLine2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
end



struct Line{T<:Number} <: AbstractLine3D{T}
    a::PVector{T}
    b::PVector{T}
end



len(line::AbstractLine) = norm(line.a - line.b)