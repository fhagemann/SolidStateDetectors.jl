# abstract type AbstractPlanarPoint{T} <: StaticArrays.FieldVector{2, T} end
# abstract type AbstractPlanarVector{T} <: StaticArrays.FieldVector{2, T} end
# 
# struct PlanarPoint{T} <: AbstractPlanarPoint{T}
#     u::T
#     v::T
# end

"""
    struct CartesianPoint{T} <: AbstractCoordinatePoint{T, Cartesian}

Describes a three-dimensional point in Cartesian coordinates.

## Fields 
* `x`: x-coordinate (in m).
* `y`: y-coordinate (in m).
* `z`: z-coordinate (in m).

See also [`CylindricalPoint`](@ref).
"""
struct CartesianPoint{T<:AbstractFloat} <: AbstractCoordinatePoint{T, Cartesian}
    x::T
    y::T
    z::T
end

#Type conversion happens here
function CartesianPoint{T}(x,y,z) where {T}
    _x = _csg_convert_args(T, x)
    _y = _csg_convert_args(T, y)
    _z = _csg_convert_args(T, z)
    CartesianPoint{T}(_x, _y, _z)
end

#Type promotion happens here
function CartesianPoint(x::TX, y::TY, z::TZ) where {TX,TY,TZ}
    eltypes = _csg_get_promoted_eltype.((TX,TY,TZ))
    T = float(promote_type(eltypes...))
    CartesianPoint{T}(x,y,z)
end

function CartesianPoint(;
    x = 0,
    y = 0,
    z = 0
)
    CartesianPoint(x,y,z)
end

function CartesianPoint{T}(;
    x = 0,
    y = 0,
    z = 0
) where {T}
    CartesianPoint{T}(x,y,z)
end

zero(PT::Type{<:AbstractCoordinatePoint{T}}) where {T} = PT(zero(T),zero(T),zero(T))

"""
    struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}

Describes a three-dimensional point in cylindrical coordinates. 

## Fields
* `r`: Radius (in m).
* `φ`: Polar angle (in rad).
* `z`: `z`-coordinate (in m).

!!! note 
    `φ == 0` corresponds to the `x`-axis in the Cartesian coordinate system.
    
See also [`CartesianPoint`](@ref).
"""
struct CylindricalPoint{T} <: AbstractCoordinatePoint{T, Cylindrical}
    r::T
    φ::T
    z::T
    CylindricalPoint{T}(r::T, φ::T, z::T) where {T} = new(r, mod(φ,T(2π)), z)
    CylindricalPoint{T}(r::Real, φ::Real, z::Real) where {T} = new(T(r), mod(T(φ),T(2π)), T(z))
end

function CylindricalPoint(pt::CartesianPoint{T})::CylindricalPoint{T} where {T}
    return CylindricalPoint{T}(hypot(pt.x, pt.y), atan(pt.y, pt.x), pt.z)
end

function CartesianPoint(pt::CylindricalPoint{T})::CartesianPoint{T} where {T}
    sφ::T, cφ::T = sincos(pt.φ)
    return CartesianPoint{T}(pt.r * cφ, pt.r * sφ, pt.z)
end

@inline CylindricalPoint(pt::CylindricalPoint) = pt
@inline CartesianPoint(pt::CartesianPoint) = pt

@inline _convert_point(pt::AbstractCoordinatePoint, ::Type{Cylindrical}) = CylindricalPoint(pt)
@inline _convert_point(pt::AbstractCoordinatePoint, ::Type{Cartesian}) = CartesianPoint(pt)

# function _Δφ(φ1::T, φ2::T)::T where {T}
#     δφ = mod(φ2 - φ1, T(2π))
#     min(δφ, T(2π) - δφ)
# end
# 
# _φNear(φ::Real, φMin::T, φMax::T) where {T} = _Δφ(T(φ),φMin) ≤ _Δφ(T(φ),φMax) ? φMin : φMax
