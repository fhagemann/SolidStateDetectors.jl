"""
    struct ImpurityScale{T, N, S, AT} <: AbstractArray{T, N}
        
Impurity scalar field of the simulation.
It is kinda an alpha map for the impurity density of semiconductors.
It can takes values between `0` and `1`:
 * `1`: The impurity density has its full value. For grid points in depleted region of the semiconductor. 
 * `]0,1[`: The impurity density is scaled down but not zero. For grid points at the edge of the depleted region which are partially depleted.
 * `0`: The impurity density is set to `0`. For grid points in undepleted regions of the semiconductor.

        
## Parametric types 
* `T`: Element type of `data`.
* `N`: Dimension of the `grid` and `data` array.  
* `S`: Coordinate system (`Cartesian` or `Cylindrical`).
* `AT`: Axes type.  
        
## Fields
* `data::Array{T, N}`: Array containing the values of the impurity scale at the discrete points of the `grid`.
* `grid::Grid{T, N, S, AT}`: [`Grid`](@ref) defining the discrete points for which the impurity scale is determined.
"""
struct ImpurityScale{T, N, S, AT} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S, AT}
end

@inline size(epot::ImpurityScale{T, N, S}) where {T, N, S} = size(epot.data)
@inline length(epot::ImpurityScale{T, N, S}) where {T, N, S} = length(epot.data)
@inline getindex(epot::ImpurityScale{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(epot.data, I...)
@inline getindex(epot::ImpurityScale{T, N, S}, i::Int) where {T, N, S} = getindex(epot.data, i)
@inline getindex(epot::ImpurityScale{T, N, S}, s::Symbol) where {T, N, S} = getindex(epot.grid, s)


function NamedTuple(epot::ImpurityScale{T, 3}) where {T}
    return (
        grid = NamedTuple(epot.grid),
        values = epot.data,
    )
end
Base.convert(T::Type{NamedTuple}, x::ImpurityScale) = T(x)
    
function ImpurityScale(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_system(grid)
    N = get_number_of_dimensions(grid)
    ImpurityScale{T, N, S, typeof(grid.axes)}(nt.values, grid)
end
Base.convert(T::Type{ImpurityScale}, x::NamedTuple) = T(x)