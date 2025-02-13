abstract type AbstractGrid{T <: SSDFloat, N} end

"""
    struct Grid{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractGrid{T, N}

Collection of `N` axes that define the dimensions of the grid needed to calculate 
[`ElectricPotential`](@ref), [`ElectricField`](@ref) or [`WeightingPotential`](@ref).

## Parametric types
* `T`: Tick type (element type) of the axes.
* `N`: Dimension of the grid.
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.

## Fields
* `axes::AT`: Tuple of length `N` containing `DiscreteAxis` for each dimension of the grid.

See also [`DiscreteAxis`](@ref).
"""
struct Grid{T, N, S <: AbstractCoordinateSystem, AT} <: AbstractGrid{T, N}
    axes::AT
end

const CartesianGrid{T, N} = Grid{T, N, Cartesian}
const CartesianGrid1D{T} = CartesianGrid{T, 1}
const CartesianGrid2D{T} = CartesianGrid{T, 2}
const CartesianGrid3D{T} = CartesianGrid{T, 3}
const CylindricalGrid{T} = Grid{T, 3, Cylindrical}

CylindricalGrid{T}(a) where {T} = Grid{T, 3, Cylindrical, typeof(a)}(a)
CartesianGrid3D{T}(a) where {T} = Grid{T, 3, Cartesian, typeof(a)}(a)

@inline Base.size(grid::Grid) = size.(grid.axes, 1)
@inline Base.size(grid::Grid, i::Int) = size(grid.axes[i], 1)
@inline Base.length(grid::Grid) = prod(size(grid))
@inline Base.getindex(grid::Grid{<:Any, N}, I::Vararg{Int, N}) where {N} = broadcast(getindex, grid.axes, I)
@inline Base.getindex(grid::Grid, i::Int) = getproperty(grid, :axes)[i]
@inline Base.getindex(grid::Grid, s::Symbol) = getindex(grid, Val{s}())

@inline Base.getproperty(grid::Grid, s::Symbol) = getproperty(grid, Val{s}())
@inline Base.getproperty(grid::Grid, ::Val{:axes}) = getfield(grid, :axes)

@inline Base.getproperty(grid::CylindricalGrid, ::Val{:axes}) = getfield(grid, :axes)
@inline Base.getproperty(grid::CylindricalGrid, ::Val{:r}) = @inbounds grid.axes[1]
@inline Base.getproperty(grid::CylindricalGrid, ::Val{:φ}) = @inbounds grid.axes[2]
@inline Base.getproperty(grid::CylindricalGrid, ::Val{:z}) = @inbounds grid.axes[3]
@inline Base.getproperty(grid::CartesianGrid3D, ::Val{:x}) = @inbounds grid.axes[1]
@inline Base.getproperty(grid::CartesianGrid3D, ::Val{:y}) = @inbounds grid.axes[2]
@inline Base.getproperty(grid::CartesianGrid3D, ::Val{:z}) = @inbounds grid.axes[3]

@inline Base.getindex(grid::CylindricalGrid, ::Val{:r}) = @inbounds grid.axes[1]
@inline Base.getindex(grid::CylindricalGrid, ::Val{:φ}) = @inbounds grid.axes[2]
@inline Base.getindex(grid::CylindricalGrid, ::Val{:z}) = @inbounds grid.axes[3]
@inline Base.getindex(grid::CartesianGrid3D, ::Val{:x}) = @inbounds grid.axes[1]
@inline Base.getindex(grid::CartesianGrid3D, ::Val{:y}) = @inbounds grid.axes[2]
@inline Base.getindex(grid::CartesianGrid3D, ::Val{:z}) = @inbounds grid.axes[3]

@inline GridPoint(grid::Grid{T, 3, Cylindrical}, inds::NTuple{3, Int}) where {T} = 
    CylindricalPoint{T}(broadcast(i -> grid.axes[i].ticks[inds[i]], (1, 2, 3)))
@inline GridPoint(grid::Grid{T, 3, Cartesian}, inds::NTuple{3, Int}) where {T} = 
    CartesianPoint{T}(broadcast(i -> grid.axes[i].ticks[inds[i]], (1, 2, 3)))

Base.sizeof(grid::Grid) = sum(sizeof.(grid.axes))

function Base.print(io::IO, grid::Grid{T, N, S}) where {T, N, S}
    print(io, "Grid{$T, $N, $S}", grid.axes)
end
function Base.println(io::IO, grid::Grid{T, N, S}) where {T, N, S}
    println(" Grid{$T, $N, $S}")
    for (i, ax) in enumerate(grid.axes)
        println(io, "  Axis $(i): ", ax)
    end
end
Base.show(io::IO, grid::Grid) = print(io, grid)
Base.show(io::IO, ::MIME"text/plain", grid::Grid) = show(io, grid)

function check_grid(grid::CylindricalGrid)::Nothing
    nr::Int, nφ::Int, nz::Int = size(grid)
    @assert iseven(nz) "GridError: Field simulation algorithm in cylindrical coordinates needs an even number of grid points in z. This is not the case. #z-ticks = $(nz)."
    @assert (iseven(nφ) || (nφ == 1)) "GridError: Field simulation algorithm in cylindrical coordinates needs an even number of grid points in φ or just one point (2D). This is not the case. #φ-ticks = $(nφ)."
    nothing
end

function check_grid(grid::CartesianGrid3D)::Nothing
    nx::Int, _, _ = size(grid)
    @assert iseven(nx) "GridError: Field simulation algorithm in cartesian coordinates needs an even number of grid points in x. This is not the case. #x-ticks = $(nx)."
    nothing
end

get_coordinate_system(grid::Grid{<:Any, <:Any, S}) where {S} = S
Base.ndims(grid::Grid{<:Any, N}) where {N} = N
get_number_of_dimensions(grid::Grid) = ndims(grid)
Base.eltype(grid::Grid{T}) where {T} = T
get_boundary_types(grid::Grid) = get_boundary_types.(grid.axes)

# Tuples with ticks to sample with differently spaced ticks
const CartesianTicksTuple{T} = NamedTuple{(:x,:y,:z), NTuple{3,Vector{T}}}
const CylindricalTicksTuple{T} = NamedTuple{(:r,:φ,:z), NTuple{3,Vector{T}}}

TicksTuple(grid::CartesianGrid3D{T}) where {T} = (x = grid.axes[1].ticks, y = grid.axes[2].ticks, z = grid.axes[3].ticks)
TicksTuple(grid::CylindricalGrid{T}) where {T} = (r = grid.axes[1].ticks, φ = grid.axes[2].ticks, z = grid.axes[3].ticks)

function Grid(nt::NamedTuple)
    if nt.coordtype == "cylindrical"
        axr::DiscreteAxis = DiscreteAxis(nt.axes.r, unit=u"m")
        axφ::DiscreteAxis = DiscreteAxis(nt.axes.phi, unit=u"rad")
        axz::DiscreteAxis = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axr.ticks[1])
        return CylindricalGrid{T}( (axr, axφ, axz) )
    elseif nt.coordtype == "cartesian"
        axx::DiscreteAxis = DiscreteAxis(nt.axes.x, unit=u"m")
        axy::DiscreteAxis = DiscreteAxis(nt.axes.y, unit=u"m")
        axz = DiscreteAxis(nt.axes.z, unit=u"m")
        T = typeof(axx.ticks[1])
        return CartesianGrid3D{T}( (axx, axy, axz) )
    else
        error("`coordtype` = $(nt.coordtype) is not valid.")
    end
end

Base.convert(T::Type{Grid}, x::NamedTuple) = T(x)

function NamedTuple(grid::CylindricalGrid{T}) where {T}
    axr::DiscreteAxis{T} = grid.axes[1]
    axφ::DiscreteAxis{T} = grid.axes[2]
    axz::DiscreteAxis{T} = grid.axes[3]
    return (
        coordtype = "cylindrical",
        ndims = 3,
        axes = (
            r   = NamedTuple(axr, unit=u"m"),
            phi = NamedTuple(axφ, unit=u"rad"),
            z   = NamedTuple(axz, unit=u"m"),
        )
    )
end
function NamedTuple(grid::CartesianGrid3D{T}) where {T}
    axx::DiscreteAxis{T} = grid.axes[1]
    axy::DiscreteAxis{T} = grid.axes[2]
    axz::DiscreteAxis{T} = grid.axes[3]
    return (
        coordtype = "cartesian",
        ndims = 3,
        axes = (
            x = NamedTuple(axx, unit=u"m"),
            y = NamedTuple(axy, unit=u"m"),
            z = NamedTuple(axz, unit=u"m"),
        )
    )
end

Base.convert(T::Type{NamedTuple}, x::Grid) = T(x)


function find_closest_gridpoint(pt::CylindricalPoint{T}, grid::CylindricalGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    return (searchsortednearest(grid.axes[1].ticks, pt.r), searchsortednearest(grid.axes[2].ticks, pt.φ), searchsortednearest(grid.axes[3].ticks, pt.z))
end
function find_closest_gridpoint(pt::CartesianPoint{T}, grid::CylindricalGrid{T})::NTuple{3,Int} where {T <: SSDFloat}
    find_closest_gridpoint(CylindricalPoint(pt),grid)
end

function find_closest_gridpoint(pt::CartesianPoint{T}, grid::CartesianGrid3D{T})::NTuple{3,Int} where {T <: SSDFloat}
    @inbounds return (searchsortednearest(grid.axes[1].ticks, pt.x), searchsortednearest(grid.axes[2].ticks, pt.y), searchsortednearest(grid.axes[3].ticks, pt.z))
end
function find_closest_gridpoint(pt::CylindricalPoint{T}, grid::CartesianGrid3D{T})::NTuple{3,Int} where {T <: SSDFloat}
    find_closest_gridpoint(CartesianPoint(pt),grid)
end

multiplicity(g::Grid) = prod(multiplicities(g))

function multiplicities(g::CylindricalGrid{T}) where {T}
    mr = one(T)
    mφ = T(2π) / width(g.axes[2].interval)
    mz = multiplicity(g.axes[3], Cartesian)
    mr, mφ, mz
end

multiplicities(g::CartesianGrid3D) = broadcast(ax -> multiplicity(ax, Cartesian), g.axes)


function voxel_widths(grid::CartesianGrid3D{T}, i1::Int, i2::Int, i3::Int) where {T} 
    wx::T = grid[1].ticks[i1 + 1] - grid[1].ticks[i1]
    wy::T = grid[2].ticks[i2 + 1] - grid[2].ticks[i2]
    wz::T = grid[3].ticks[i3 + 1] - grid[3].ticks[i3]
    wx, wy, wz
end
function voxel_widths(grid::CylindricalGrid{T}, i1::Int, i2::Int, i3::Int) where {T} 
    wr::T =  grid[1].ticks[i1 + 1] - grid[1].ticks[i1]
    wφ::T = (grid[2].ticks[i2 + 1] - grid[2].ticks[i2]) * (grid[1].ticks[i1 + 1] + grid[1].ticks[i1])/2
    wz::T =  grid[3].ticks[i3 + 1] - grid[3].ticks[i3]
    wr, wφ, wz
end

voxel_volume(grid::CylindricalGrid{T}, i1::Int, i2::Int, i3::Int, w1::T, w2::T, w3::T) where {T} = 
    (grid[2].ticks[i2 + 1] - grid[2].ticks[i2]) * w3 * (grid[1].ticks[i1 + 1]^2 - grid[1].ticks[i1]^2) / 2  
voxel_volume(grid::CartesianGrid3D{T}, i1::Int, i2::Int, i3::Int, w1::T, w2::T, w3::T) where {T} =
    w1 * w2 * w3

    
function get_extended_midpoints_grid(grid::Grid{T,3}) where {T}
    ticks = broadcast(i -> midpoints(get_extended_ticks(grid.axes[i])), (1,2,3))
    axes = broadcast(i -> typeof(grid.axes[i])(grid.axes[i].interval, ticks[i]) , (1,2,3))
    typeof(grid)(axes)
end

in(pt::CartesianPoint{T}, grid::Grid{T, 3, Cartesian}) where {T} =
    all(broadcast(i -> pt[i] in grid.axes[i].interval, (1, 2, 3)))
in(pt::CylindricalPoint{T}, grid::Grid{T, 3, Cylindrical}) where {T} =
    all(broadcast(i -> pt[i] in grid.axes[i].interval, (1, 2, 3)))
in(pt::CartesianPoint{T}, grid::Grid{T, 3, Cylindrical}) where {T} =
    in(CylindricalPoint(pt), grid)