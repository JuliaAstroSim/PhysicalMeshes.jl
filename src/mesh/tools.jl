"""
$(TYPEDSIGNATURES)
Enlarge the array by duplicating elements in situ, or by interpolation.

Often used interpolation algorithms: `LinearInterpolation`, `BSplineInterpolation`, `CubicSplineInterpolation`

# Examples
```
julia> zoom([1,2],2)
4-element Vector{Int64}:
 1
 1
 2
 2

julia> zoom([1 2; 3 4],2)
4×4 Matrix{Int64}:
 1  1  2  2
 1  1  2  2
 3  3  4  4
 3  3  4  4

julia> zoom(reshape([1],1,1,1),2)
2×2×2 Array{Int64, 3}:
[:, :, 1] =
 1  1
 1  1
  
[:, :, 2] =
 1  1
 1  1

julia> using Interpolations

julia> zoom([1,5,3],2, interp = LinearInterpolation)
 6-element Vector{Float64}:
  1.0
  2.5999999999999996
  4.200000000000001
  4.6
  3.8
  3.0
```
"""
function zoom(a::AbstractArray{T,1}, ratio::Int; interp = nothing) where T
    Nx = size(a, 1)
    if isnothing(interp)
        out = Array{T,1}(undef, Nx*ratio)
        for i in 1:Nx*ratio
            @inbounds out[i] = a[div(i-1,ratio)+1]
        end
    else
        out = zeros(Nx*ratio)
        xs = 1:Nx
        itp = interp(xs, a)

        for i in 1:Nx*ratio
            x = (Nx-1)/(ratio*Nx-1)*i + (ratio-1)*Nx/(ratio*Nx-1)
            x = x > Nx ? Nx : x
            @inbounds out[i] = itp(x)
        end
    end
    out
end

"""
$(TYPEDSIGNATURES)
"""
function zoom(a::AbstractArray{T,2}, ratio::Int; interp = nothing) where T
    Nx, Ny = size(a)
    if isnothing(interp)
        out = Array{T,2}(undef, Nx*ratio, Ny*ratio)
        for j in 1:Ny*ratio
            for i in 1:Nx*ratio
                @inbounds out[i, j] = a[div(i-1,ratio)+1, div(j-1,ratio)+1]
            end
        end
    else
        out = zeros(Nx*ratio, Ny*ratio)
        xs = 1:Nx
        ys = 1:Ny
        itp = interp((xs,ys), a)

        for j in 1:Ny*ratio
            for i in 1:Nx*ratio
                x = (Nx-1)/(ratio*Nx-1)*i + (ratio-1)*Nx/(ratio*Nx-1)
                x = x > Nx ? Nx : x
                y = (Ny-1)/(ratio*Ny-1)*j + (ratio-1)*Ny/(ratio*Ny-1)
                y = y > Ny ? Ny : y
                @inbounds out[i, j] = itp(x, y)
            end
        end
    end
    out
end

"""
$(TYPEDSIGNATURES)
"""
function zoom(a::AbstractArray{T,3}, ratio::Int; interp = nothing) where T
    Nx, Ny, Nz = size(a)
    if isnothing(interp)
        out = Array{T,3}(undef, Nx*ratio, Ny*ratio, Nz*ratio)
        for k in 1:Nz*ratio
            for j in 1:Ny*ratio
                for i in 1:Nx*ratio
                    @inbounds out[i, j, k] = a[div(i-1,ratio)+1, div(j-1,ratio)+1, div(k-1,ratio)+1]
                end
            end
        end
    else
        out = zeros(Nx*ratio, Ny*ratio, Nz*ratio)
        xs = 1:Nx
        ys = 1:Ny
        zs = 1:Nz
        itp = interp((xs,ys,zs), a)

        for k in 1:Nz*ratio
            for j in 1:Ny*ratio
                for i in 1:Nx*ratio
                    x = (Nx-1)/(ratio*Nx-1)*i + (ratio-1)*Nx/(ratio*Nx-1)
                    x = x > Nx ? Nx : x
                    y = (Ny-1)/(ratio*Ny-1)*j + (ratio-1)*Ny/(ratio*Ny-1)
                    y = y > Ny ? Ny : y
                    z = (Nz-1)/(ratio*Nz-1)*k + (ratio-1)*Nz/(ratio*Nz-1)
                    z = z > Nz ? Nz : z
                    @inbounds out[i, j, k] = itp(x, y, z)
                end
            end
        end
    end
    out
end

"""
$(TYPEDSIGNATURES)
shrink the array size by simply selecting elements whose index can divide the ratio exactly
# Examples
```
julia> shrink([1,2,3,4,5,6],2)
3-element Vector{Int64}:
 2
 4
 6

julia> shrink([1,2,3,4,5,6],3)
2-element Vector{Int64}:
 3
 6

julia> shrink(ones(4,4),2)
2×2 Matrix{Float64}:
 1.0  1.0
 1.0  1.0

julia> shrink(ones(6,6,6),3)
2×2×2 Array{Float64, 3}:
[:, :, 1] =
 1.0  1.0
 1.0  1.0

[:, :, 2] =
 1.0  1.0
 1.0  1.0
```
"""
function shrink(a::AbstractArray{T,1}, ratio::Int) where T
    return a[iszero.(mod.(eachindex(a), ratio))]
end

"""
$(TYPEDSIGNATURES)
"""
function shrink(a::AbstractArray{T,2}, ratio::Int) where T
    Nx, Ny = size(a)
    out = Array{T,2}(undef, div(Nx, ratio), div(Ny, ratio))
    for j in 1:div(Ny,ratio)
        for i in 1:div(Nx,ratio)
            @inbounds out[i, j] = a[i*ratio, j*ratio]
        end
    end
    out
end

"""
$(TYPEDSIGNATURES)
"""
function shrink(a::AbstractArray{T,3}, ratio::Int) where T
    Nx, Ny, Nz = size(a)
    out = Array{T,3}(undef, div(Nx,ratio), div(Ny,ratio), div(Nz,ratio))
    for k in 1:div(Nz,ratio)
        for j in 1:div(Ny,ratio)
            for i in 1:div(Nx,ratio)
                @inbounds out[i, j, k] = a[i*ratio, j*ratio, k*ratio]
            end
        end
    end
    out
end

function mass_center(m::MeshCartesianStatic)
    rho_total = sum(m.rho)
    if iszero(rho_total)
        return mean(m.pos)
    else
        return mean(m.rho .* m.pos) / rho_total
    end
end

function total_mass(m::MeshCartesianStatic)
    return sum(m.rho .* prod(m.config.Δ))
end