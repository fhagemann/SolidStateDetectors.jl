abstract type AbstractChargeCloud end

struct PointCharge{T <: SSDFloat, S} <: AbstractChargeCloud 
    charge::T
    pos::AbstractCoordinatePoint{T, 3, S}
end

struct NBodyChargeCloud{T <: SSDFloat, S} <: AbstractChargeCloud 
    # To be done...
end

struct Tetrahedron{T}
    points::Vector{CartesianPoint{T}}
end

function Tetrahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[center + CartesianVector{T}(0,0,length)]
    for φ in (0,120,240)
        push!(points, center + length * CartesianVector{T}(sqrt(8)/3*cosd(φ), sqrt(8)/3*sind(φ), -1/3))
    end
    Tetrahedron{T}( points )
end

struct Hexahedron{T}
    points::Vector{CartesianPoint{T}}
end

function Hexahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[]
    for x in (-1,1)
        for y in (-1,1)
            for z in (-1,1)
                push!(points, center + length * sqrt(T(1/3)) * CartesianVector{T}(x,y,z))
            end
        end
    end
    Hexahedron{T}( points )
end

struct Octahedron{T}
    points::Vector{CartesianPoint{T}}
end

function Octahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[
        center + length * CartesianVector{T}(1,0,0),
        center + length * CartesianVector{T}(-1,0,0),
        center + length * CartesianVector{T}(0,1,0),
        center + length * CartesianVector{T}(0,-1,0),
        center + length * CartesianVector{T}(0,0,1),
        center + length * CartesianVector{T}(0,0,-1)
    ]
    Octahedron{T}( points )
end

struct Dodecahedron{T}
    points::Vector{CartesianPoint{T}}
end

function Dodecahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[
        center + sqrt(1/3) * length * CartesianVector{T}(-1,-1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,-1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,-1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,-1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,1,1),
    ]
    for g in (-Base.MathConstants.golden, Base.MathConstants.golden)
        for invg in (-inv(Base.MathConstants.golden), inv(Base.MathConstants.golden))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(invg, g, 0))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(0, invg, g))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(g, 0, invg))
        end
    end
    Dodecahedron{T}( points )
end

struct Icosahedron{T}
    points::Vector{CartesianPoint{T}}
end

function Icosahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[]
    for i in (-1,1)
        for g in (-Base.MathConstants.golden, Base.MathConstants.golden)
            push!(points, center + length * CartesianVector{T}(i,0,g) / sqrt(T(2 + Base.MathConstants.golden)))
            push!(points, center + length * CartesianVector{T}(g,i,0) / sqrt(T(2 + Base.MathConstants.golden)))
            push!(points, center + length * CartesianVector{T}(0,g,i) / sqrt(T(2 + Base.MathConstants.golden)))
        end
    end
    Icosahedron{T}( points )
end

include("plot_recipes.jl")
