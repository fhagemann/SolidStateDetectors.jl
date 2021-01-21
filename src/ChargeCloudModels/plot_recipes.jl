@recipe function f(t::PointCharge{T}) where {T}

    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "PointCharge"
        t.points
    end
end


@recipe function f(t::Tetrahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Tetrahedron"
        t.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(1,2,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(3,1,4)]
        end
    end
end


@recipe function f(i::Hexahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Hexahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,6,5,1,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(4,8,7,3)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,7)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,8)]
        end
    end
end


@recipe function f(i::Octahedron{T}) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Octahedron"
        i.points
    end
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :line3d
        label --> ""
        linewidth --> 1
        i.points[vcat(1,3,2,4,1,5,2,6,3,5,4,6,1)]
    end
end


@recipe function f(i::Dodecahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Dodecahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,9,2,14,5,15,3,11,1,10,4,12,6,16,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,19,16)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(19,8,18,15)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,13,10)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(4,17,7,13)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(7,18)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(8,20,17)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,20)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(9,12)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(11,14)]
        end
    end
end


@recipe function f(i::Icosahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Icosahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(7,6,1,2,3,5,7,3,1,7,11,5,9,2,8,1)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,9,10,4,2,9,4,8,6,11,10,12,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,10,12,8)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,12,11)]
        end
    end
end


get_vertices(::Type{PointCharge{T}}) where {T} = 1
get_vertices(::Type{Tetrahedron{T}}) where {T} = 4
get_vertices(::Type{Hexahedron{T}}) where {T} = 8
get_vertices(::Type{Octahedron{T}}) where {T} = 6
get_vertices(::Type{Dodecahedron{T}}) where {T} = 20
get_vertices(::Type{Icosahedron{T}}) where {T} = 12


@recipe function f(nbc::NBodyChargeCloud{T}, connect = true) where {T}
    
    seriescolor --> :blue
    markersize --> 10
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> markersize
        markerstrokewidth --> 0
        label --> "NBodyChargeCloud"
        []
    end
    
    points = nbc.points
    vertex_no = 0
    for (shell, shell_structure) in enumerate(nbc.shell_structure)
        vertices = get_vertices(shell_structure)
        @series begin
            markersize --> markersize * exp(-(shell - 1))
            label --> ""
            seriescolor --> seriescolor
            shell_structure(points[vertex_no+1:vertex_no+vertices])
        end
        vertex_no += vertices
    end
end