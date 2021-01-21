abstract type AbstractChargeCloud end

struct NBodyChargeCloud{T} <: AbstractChargeCloud
    points::Vector{CartesianPoint{T}}
    charges::Vector{T}
    shell_structure::Vector{Type{<:AbstractChargeCloud}}
end

struct PointCharge{T} <: AbstractChargeCloud
    points::SVector{1, CartesianPoint{T}}
end

function PointCharge(center::Vector{CartesianPoint{T}}, length::T = T(0)) where {T} 
    PointCharge{T}(center)
end

struct Tetrahedron{T} <: AbstractChargeCloud
    points::SVector{4, CartesianPoint{T}}
end

function Tetrahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[center + CartesianVector{T}(0,0,length)]
    for φ in (0,120,240)
        push!(points, center + length * CartesianVector{T}(sqrt(8)/3*cosd(φ), sqrt(8)/3*sind(φ), -1/3))
    end
    Tetrahedron{T}( points )
end

struct Hexahedron{T} <: AbstractChargeCloud
    points::SVector{8, CartesianPoint{T}}
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

struct Octahedron{T} <: AbstractChargeCloud
    points::SVector{6, CartesianPoint{T}}
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

struct Dodecahedron{T} <: AbstractChargeCloud
    points::SVector{20, CartesianPoint{T}}
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

struct Icosahedron{T} <: AbstractChargeCloud
    points::SVector{12, CartesianPoint{T}}
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
include("ParticleTypes.jl")


function create_charge_cloud(center::CartesianPoint{T}, charge::T, particle_type::Type{PT} = Gamma;
        radius::T = radius_guess(charge, particle_type), number_of_shells::Int = 2, shell_structure = Dodecahedron
    )::NBodyChargeCloud{T} where {T, PT <: ParticleType}
    
    points::Vector{CartesianPoint{T}} = CartesianPoint{T}[center]
    charges::Vector{T} = T[charge]
    shell_structures::Vector{Type} = [PointCharge{T}]
    
    n_shell = 1
    while n_shell <= number_of_shells
        shell = shell_structure(center, n_shell * radius).points
        points = vcat(points, shell)
        charges = vcat(charges, [exp(-n_shell^2 / 2) for i in 1:length(shell)])
        push!(shell_structures, shell_structure{T})
        n_shell += 1
    end
    return NBodyChargeCloud{T}( points, charges./sum(charges) * charge, shell_structures )
end


