function LineSegments(t::SolidStateDetectors.Tube{T}; n::Int = 30)::Vector{SolidStateDetectors.LineSegment{T, 3, :cartesian}} where {T <: SolidStateDetectors.SSDFloat}
    ls = SolidStateDetectors.LineSegment{T, 3, :cartesian}[]
    φs = range(t.φ_interval.left, length = n+1, stop = t.φ_interval.right)
    for r in (t.r_interval.left == 0 ? [t.r_interval.right] : [t.r_interval.left, t.r_interval.right])
        for z in [t.z_interval.left, t.z_interval.right]
            pts = [CartesianPoint{T}(r * sin(φ), r * cos(φ), z) for φ in φs]
            for i in 1:n
                push!(ls, SolidStateDetectors.LineSegment(pts[i], pts[i+1]))
            end
        end
    end
    for r in [t.r_interval.left, t.r_interval.right]
        if r != 0
            for φ in ((φs[end] - φs[1] ≈ 2π) ? [φs[1]] : [φs[1], φs[end]])
                push!(ls, SolidStateDetectors.LineSegment(
                    CartesianPoint{T}(r * sin(φ), r * cos(φ), t.z_interval.left),
                    CartesianPoint{T}(r * sin(φ), r * cos(φ), t.z_interval.right)))
            end
        end
    end
    for φ in ((φs[end] - φs[1] ≈ 2π) ? [] : [φs[1], φs[end]])
        for z in [t.z_interval.left, t.z_interval.right]
            push!(ls, SolidStateDetectors.LineSegment(
                CartesianPoint{T}(t.r_interval.left * sin(φ),  t.r_interval.left * cos(φ),  z),
                CartesianPoint{T}(t.r_interval.right * sin(φ), t.r_interval.right * cos(φ), z)))
        end
    end
    return ls
end

@recipe function f(t::Tube{T}; n = 30, seriescolor = :green) where {T}
    linewidth --> 2
    @series begin
        seriescolor --> seriescolor
        label --> "Tube"
        []
    end
    label := ""
    seriescolor := seriescolor
    LineSegments(t, n = n)
end



function LineSegments(c::SolidStateDetectors.Cone{T}; n::Int = 30)::Vector{SolidStateDetectors.LineSegment{T, 3, :cartesian}} where {T <: SolidStateDetectors.SSDFloat}
    ls = SolidStateDetectors.LineSegment{T, 3, :cartesian}[]
    φs = range(c.φStart, length = n+1, stop = c.φStop)
    for r in [c.rStart1, c.rStop1]
        for z in [c.zStart]
            pts = [CartesianPoint{T}(r * sin(φ), r * cos(φ), z) for φ in φs]
            for i in 1:n
                push!(ls, SolidStateDetectors.LineSegment(pts[i], pts[i+1]))
            end
        end
    end
    for r in [c.rStart2, c.rStop2]
        for z in [c.zStop]
            pts = [CartesianPoint{T}(r * sin(φ), r * cos(φ), z) for φ in φs]
            for i in 1:n
                push!(ls, SolidStateDetectors.LineSegment(pts[i], pts[i+1]))
            end
        end
    end
    for φ in ((φs[end] - φs[1] ≈ 2π) ? [φs[1]] : [φs[1], φs[end]])
        push!(ls, SolidStateDetectors.LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ), c.rStart1 * cos(φ), c.zStart),
            CartesianPoint{T}(c.rStart2 * sin(φ), c.rStart2 * cos(φ), c.zStop)))
        push!(ls, SolidStateDetectors.LineSegment(
            CartesianPoint{T}(c.rStop1 * sin(φ), c.rStop1 * cos(φ), c.zStart),
            CartesianPoint{T}(c.rStop2 * sin(φ), c.rStop2 * cos(φ), c.zStop)))
    end
    for φ in ((φs[end] - φs[1] ≈ 2π) ? [φs[1]] : [φs[1], φs[end]])
        push!(ls, SolidStateDetectors.LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ),  c.rStart1 * cos(φ),  c.zStart),
            CartesianPoint{T}(c.rStop1 * sin(φ),   c.rStop1 * cos(φ), c.zStart)))
        push!(ls, SolidStateDetectors.LineSegment(
            CartesianPoint{T}(c.rStart2 * sin(φ),  c.rStart2 * cos(φ),  c.zStop),
            CartesianPoint{T}(c.rStop2 * sin(φ),   c.rStop2 * cos(φ), c.zStop)))
    end
    return ls
end

@recipe function f(c::Cone{T}; n = 30, seriescolor = :orange) where {T}
    linewidth --> 2
    @series begin
        seriescolor --> seriescolor
        label --> "Cone"
        []
    end
    seriescolor := seriescolor
    label := ""
    LineSegments(c, n = n)
end


#=
@recipe function f(Vol::Tube{T}, coloring = missing) where T <: SSDFloat
    rStart = Vol.r_interval.left
    rStop = Vol.r_interval.right
    φStart = Vol.φ_interval.left
    φStop = Vol.φ_interval.right
    zStart = Vol.z_interval.left
    zStop = Vol.z_interval.right

    ismissing(Vol.translate) ? translate = [0.0,0.0,0.0] : translate = Vol.translate
    ismissing(coloring) ? c --> :green : c --> coloring
    @series begin
        label --> "Tube"
        partialcircle_3d(rStop,φStart,φStop,[0,0,zStart] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStop,φStart,φStop,[0,0,zStop] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStart,φStart,φStop,[0,0,zStart] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStart,φStart,φStop,[0,0,zStop] .+ translate)
    end
    ## Vertical Lines

    if !iszero(rStart)
        @series begin
            label:= ""
            line_3d(rStart, rStart, φStart, φStart, zStart, zStop, translate = translate)
        end
    end
    if !iszero(rStart)
        @series begin
            label:= ""
            line_3d(rStart, rStart, φStop, φStop, zStart, zStop, translate = translate)
        end
    end
    @series begin
        label:= ""
        line_3d(rStop, rStop, φStart, φStart, zStart, zStop, translate = translate)
    end
    @series begin
        label:= ""
        line_3d(rStop, rStop, φStop, φStop, zStart, zStop, translate = translate)
    end

    ##Horizontal Lines

    if !isapprox((φStop - φStart)%2π , 0.0,atol=0.00001)
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStart, φStart, zStart, zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStop, φStop, zStart, zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStart, φStart, zStop, zStop, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStop, φStop, zStop, zStop, translate = translate)
        end
    end
end


@recipe function f(vol::Cone{T},n_aux_lines = 0, coloring = missing) where T
    ismissing(vol.translate) ? translate = [0.0,0.0,0.0] : translate = vol.translate
    ismissing(coloring) ? c --> :orange : c --> coloring
    @series begin
        label --> "Cone"
        partialcircle_3d(vol.rStop1,vol.φStart,vol.φStop,[0,0,vol.zStart]+translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(vol.rStop2,vol.φStart,vol.φStop,[0,0,vol.zStop]+translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(vol.rStart1,vol.φStart,vol.φStop,[0,0,vol.zStart]+translate)
    end
    @series begin
        label:= ""
        c-->coloring
        partialcircle_3d(vol.rStart2,vol.φStart,vol.φStop,[0,0,vol.zStop]+translate)
    end

    ## Vertical Lines


    @series begin
        label:= ""
        line_3d(vol.rStart1,vol.rStart2,vol.φStart,vol.φStart,vol.zStart,vol.zStop, translate = translate)
    end

    @series begin
        label:= ""
        line_3d(vol.rStart1,vol.rStart2,vol.φStop,vol.φStop,vol.zStart,vol.zStop, translate = translate)
    end


    @series begin
        label:= ""
        line_3d(vol.rStop1,vol.rStop2,vol.φStart,vol.φStart,vol.zStart,vol.zStop, translate = translate)
    end
    @series begin
        label:= ""
        line_3d(vol.rStop1,vol.rStop2,vol.φStop,vol.φStop,vol.zStart,vol.zStop, translate = translate)
    end

    ##Horizontal Lines

    if !isapprox((vol.φStop - vol.φStart)%2π , 0.0,atol=0.00001)
        @series begin
            label:= ""
            line_3d(vol.rStart1,vol.rStop1,vol.φStart,vol.φStart,vol.zStart,vol.zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart1,vol.rStop1,vol.φStop,vol.φStop,vol.zStart,vol.zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart2,vol.rStop2,vol.φStart,vol.φStart,vol.zStop,vol.zStop, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart2,vol.rStop2,vol.φStop,vol.φStop,vol.zStop,vol.zStop, translate = translate)
        end
    end

end
function partialcircle_3d(radius,phiStart,phiStop,Translate::AbstractVector;nSteps=400)
    phirange = mylinspace(phiStart,phiStop,nSteps)

    x::Vector{AbstractFloat}=map(x->radius*cos.(x) , phirange)
    y::Vector{AbstractFloat}=map(x->radius*sin.(x) , phirange)
    z::Vector{AbstractFloat}=map(x->Translate[3], phirange)
    x.+=Translate[1]
    y.+=Translate[2]

    return x,y,z
end
=#
