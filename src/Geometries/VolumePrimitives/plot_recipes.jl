function LineSegments(t::Tube{T})::Vector{AbstractLine{T,3,:cartesian}} where {T <: SSDFloat}
    ls = AbstractLine{T, 3, :cartesian}[]
    translate::CartesianVector{T} = ismissing(t.translate) ? CartesianVector{T}([0, 0, 0]) : t.translate
    for r in (t.r_interval.left == 0 ? [t.r_interval.right] : [t.r_interval.left, t.r_interval.right])
        for z in [t.z_interval.left, t.z_interval.right]
            push!(ls, PartialCircle(r, t.φ_interval.left, t.φ_interval.right, translate + CartesianVector{T}([0, 0, z])))
        end
    end
    for r in [t.r_interval.left, t.r_interval.right]
        if r != 0
            for φ in ((t.φ_interval.right - t.φ_interval.left ≈ 2π) ? [t.φ_interval.left] : [t.φ_interval.left, t.φ_interval.right])
                push!(ls, LineSegment(
                    CartesianPoint{T}(r * cos(φ), r * sin(φ), t.z_interval.left) + translate,
                    CartesianPoint{T}(r * cos(φ), r * sin(φ), t.z_interval.right) + translate))
            end
        end
    end
    for φ in ((t.φ_interval.right - t.φ_interval.left ≈ 2π) ? [] : [t.φ_interval.left, t.φ_interval.right])
        for z in [t.z_interval.left, t.z_interval.right]
            push!(ls, LineSegment(
                CartesianPoint{T}(t.r_interval.left * cos(φ), t.r_interval.left * sin(φ), z) + translate,
                CartesianPoint{T}(t.r_interval.right * cos(φ), t.r_interval.right * sin(φ), z) + translate))
        end
    end
    return ls
end

@recipe function f(t::Tube{T}; n = 30, seriescolor = :green) where {T}
    linewidth --> 2
    n --> n
    @series begin
        seriescolor --> seriescolor
        label --> "Tube"
        []
    end
    label := ""
    seriescolor := seriescolor
    LineSegments(t)
end


function LineSegments(c::Cone{T})::Vector{AbstractLine{T, 3, :cartesian}} where {T <: SSDFloat}
    ls = AbstractLine{T, 3, :cartesian}[]
    translate::CartesianVector{T} = ismissing(c.translate) ? CartesianVector{T}([0, 0, 0]) : c.translate
    for r in [c.rStart1, c.rStop1]
        push!(ls, PartialCircle(r, c.φStart, c.φStop, translate + CartesianVector{T}([0, 0, c.zStart]))) 
    end
    for r in [c.rStart2, c.rStop2]
        push!(ls, PartialCircle(r, c.φStart, c.φStop, translate + CartesianVector{T}([0, 0, c.zStop])))

    end
    for φ in ((c.φStop - c.φStart ≈ 2π) ? [c.φStart] : [c.φStart, c.φStop])
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ), c.rStart1 * cos(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStart2 * sin(φ), c.rStart2 * cos(φ), c.zStop) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStop1 * sin(φ), c.rStop1 * cos(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStop2 * sin(φ), c.rStop2 * cos(φ), c.zStop) + translate))
    end
    for φ in ((c.φStop - c.φStart ≈ 2π) ? [c.φStart] : [c.φStart, c.φStop])
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ),  c.rStart1 * cos(φ),  c.zStart) + translate,
            CartesianPoint{T}(c.rStop1 * sin(φ),   c.rStop1 * cos(φ), c.zStart) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart2 * sin(φ),  c.rStart2 * cos(φ),  c.zStop) + translate,
            CartesianPoint{T}(c.rStop2 * sin(φ),   c.rStop2 * cos(φ), c.zStop) + translate))
    end
    return ls
end

@recipe function f(c::Cone{T}; n = 30, seriescolor = :orange) where {T}
    linewidth --> 2
    n --> n
    @series begin
        seriescolor --> seriescolor
        label --> "Cone"
        []
    end
    seriescolor := seriescolor
    label := ""
    LineSegments(c)
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
