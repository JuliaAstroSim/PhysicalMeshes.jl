struct Line2D{T<:Union{Number, Quantity}} <: AbstractLine2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
end



struct Line{T<:Union{Number, Quantity}} <: AbstractLine3D{T}
    a::PVector{T}
    b::PVector{T}
end



len(line::AbstractLine) = norm(line.a - line.b)