struct DriftPath{T <: SSDFloat}
    path::Vector{<:AbstractCoordinatePoint{RealQuantity}}
    timestamps::Vector{RealQuantity}
end

struct EHDriftPath{T <: SSDFloat, TT <: RealQuantity}
    e_path::Vector{<:AbstractCoordinatePoint{T}}
    h_path::Vector{<:AbstractCoordinatePoint{T}}
    timestamps_e::Vector{TT}
    timestamps_h::Vector{TT}
end

_common_length(dp::EHDriftPath{T} where {T <: SSDFloat})::Int =
    max(length(dp.timestamps_e), length(dp.timestamps_h))
_common_length(dps::Vector{EHDriftPath{T}} where {T <: SSDFloat})::Int =
    maximum(_common_length.(dps))

function _common_time(dp::EHDriftPath{T, TT})::TT where {T <: SSDFloat, TT<:RealQuantity}
    max(last(dp.timestamps_e), last(dp.timestamps_h))
end
_common_time(dps::Vector{<:EHDriftPath}) =
maximum(_common_time.(dps))

function _common_timestamps(dp::Union{<:EHDriftPath{T}, Vector{<:EHDriftPath{T}}}, Δt) where {T}
    range(zero(Δt), step = Δt, stop = typeof(Δt)(_common_time(dp)) + Δt)
end

function get_velocity_vector(interpolation_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, point::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}
    return CartesianVector{T}(interpolation_field(point.x, point.y, point.z))
end

@inline function get_velocity_vector(interpolated_vectorfield, point::CylindricalPoint{T}) where {T <: SSDFloat}
    return CartesianVector{T}(interpolated_vectorfield(point.r, point.φ, point.z))
end


function _drift_charges(detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::Vector{CartesianPoint{T}}, energies::Vector{T},
                        electric_field::Interpolations.Extrapolation{<:SVector{3}, 3},
                        cdm::AbstractChargeDriftModel{T},
                        Δt::RQ, diffusion::Bool = false, self_repulsion::Bool = false; 
                        max_nsteps::Int = 2000, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat, RQ <: RealQuantity}

    drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(starting_points))
    n_events::Int = length(starting_points)
    dt::T = T(to_internal_units(internal_time_unit, Δt))

    drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_events, max_nsteps)
    drift_path_h::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_events, max_nsteps)
    timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)
    timestamps_h::Vector{T} = Vector{T}(undef, max_nsteps)
    n_e::Int = _drift_charge!(drift_path_e, timestamps_e, detector, point_types, grid, starting_points, energies./T(ustrip(detector.semiconductors[1].material.E_ionisation)), dt, electric_field, cdm, Electron, diffusion, self_repulsion, verbose = verbose)
    n_h::Int = _drift_charge!(drift_path_h, timestamps_h, detector, point_types, grid, starting_points, energies./T(ustrip(detector.semiconductors[1].material.E_ionisation)), dt, electric_field, cdm, Hole, diffusion, self_repulsion, verbose = verbose)
    
    for i in eachindex(starting_points)
        drift_paths[i] = EHDriftPath{T, T}( drift_path_e[i,1:n_e], drift_path_h[i,1:n_h], timestamps_e[1:n_e], timestamps_h[1:n_h] )
    end

    return drift_paths
end

#=
function _drift_charge( detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                       starting_point::CartesianPoint{<:SSDFloat},
                       velocity_field_e::Interpolations.Extrapolation{SVector{3, T}, 3},
                       velocity_field_h::Interpolations.Extrapolation{SVector{3, T}, 3},
                       cdm::AbstractChargeDriftModel{T},
                       Δt::RealQuantity, diffusion::Bool = false, self_repulsion::Bool = false;
                       max_nsteps::Int = 2000, verbose::Bool = true)::Vector{EHDriftPath{T}} where {T <: SSDFloat}
    return _drift_charges(detector, grid, CartesianPoint{T}.(point_types), [starting_point], velocity_field_e, velocity_field_h, cdm, T(Δt.val) * unit(Δt), diffusion, self_repulsion, max_nsteps = max_nsteps, verbose = verbose)
end
=#

@inline _convert_vector(pt::CartesianPoint, ::Val{:cylindrical}) = CylindricalPoint(pt)
@inline _convert_vector(pt::CartesianPoint, ::Val{:cartesian}) = pt

function modulate_surface_drift(p::CartesianVector{T})::CartesianVector{T} where {T <: SSDFloat}
    return p
end

function modulate_driftvector(sv::CartesianVector{T}, cp::CartesianPoint{T}, vdv::Vector{AbstractVirtualVolume{T}})::CartesianVector{T} where {T <: SSDFloat}
    for i in eachindex(vdv)
        if in(cp, vdv[i])
            return modulate_driftvector(sv, cp, vdv[i])
        end
    end
    return sv
end

@inline function _is_next_point_in_det(pt_car::CartesianPoint{T}, pt_cyl::CylindricalPoint{T}, det::SolidStateDetector{T, :cylindrical}, point_types::PointTypes{T, 3, :cylindrical})::Bool where {T <: SSDFloat}
    pt_cyl in point_types || pt_cyl in det
end
@inline function _is_next_point_in_det(pt_car::CartesianPoint{T}, pt_cyl::CylindricalPoint{T}, det::SolidStateDetector{T, :cartesian}, point_types::PointTypes{T, 3, :cartesian})::Bool where {T <: SSDFloat}
    pt_car in point_types || pt_car in det
end


is_zero_vector(v::CartesianVector{T}) where {T <: SSDFloat} = (v == CartesianVector{T}(0,0,0))

function _set_to_zero_vector!(v::Vector{CartesianVector{T}})::Nothing where {T <: SSDFloat}
    for n in eachindex(v)
        v[n] = CartesianVector{T}(0,0,0)
    end
    nothing
end

function _add_fieldvector_drift!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, 
        done::Vector{Bool}, velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, det::SolidStateDetector{T, S})::Nothing where {T <: SSDFloat, S}
    for n in eachindex(step_vectors)
        if !done[n]
            step_vectors[n] += get_velocity_vector(velocity_field, _convert_vector(current_pos[n], Val(S)))
            if is_zero_vector(geom_round.(step_vectors[n]))
                done[n] = true
                #current_pos[n] += step_vectors[n]
                #drift_path[n,istep] = current_pos[n]
            end
        end
    end
    nothing
end

function _add_fieldvector_diffusion!(step_vectors::Vector{CartesianVector{T}}, length::T = T(0.5e3))::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        sinθ::T, cosθ::T = sincos(T(rand())*T(2π))
        sinφ::T, cosφ::T = sincos(T(rand())*T(π))
        step_vectors[n] += CartesianVector{T}( length * cosφ * sinθ, length * sinφ * sinθ, length * cosθ )
    end
    nothing 
end

abstract type Electron end
abstract type Hole end

function _add_fieldvector_selfrepulsion!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, charges::Vector{T}, ϵr::T, ::Type{Electron})::Nothing where {T <: SSDFloat}
    #TO DO: ignore charges that are already collected (not trapped though!)
    for n in eachindex(step_vectors)
        if done[n] continue end
        for m in eachindex(step_vectors)
            if done[m] continue end
            if m > n
                direction::CartesianVector{T} = current_pos[n] .- current_pos[m]
                tmp::T = elementary_charge * inv(4 * pi * ϵ0 * ϵr * sum(direction.^2))
                step_vectors[n] -= charges[m] * tmp * normalize(direction)
                step_vectors[m] += charges[n] * tmp * normalize(direction)
            end
        end
    end
    nothing
end


function _add_fieldvector_selfrepulsion!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, done::Vector{Bool}, charges::Vector{T}, ϵr::T, ::Type{Hole})::Nothing where {T <: SSDFloat}
    #TO DO: ignore charges that are already collected (not trapped though!)
    for n in eachindex(step_vectors)
        if done[n] continue end
        for m in eachindex(step_vectors)
            if done[m] continue end
            if m > n
                direction::CartesianVector{T} = current_pos[n] .- current_pos[m]
                tmp::T = elementary_charge * inv(4 * pi * ϵ0 * ϵr * sum(direction.^2))
                step_vectors[n] += charges[m] * tmp * normalize(direction)
                step_vectors[m] -= charges[n] * tmp * normalize(direction)
            end
        end
    end
    nothing
end


function _modulate_fieldvector!(step_vectors::Vector{CartesianVector{T}}, current_pos::Vector{CartesianPoint{T}}, vdv::Vector{V})::Nothing where {T <: SSDFloat, V <: AbstractVirtualVolume{T}}
    for n in eachindex(step_vectors)
        step_vectors[n] = modulate_driftvector(step_vectors[n], current_pos[n], vdv)
    end
end

function _get_step_vectors!(step_vectors::Vector{CartesianVector{T}}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Electron})::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        step_vectors[n] = getVe(SVector{3,T}(step_vectors[n]), cdm) * Δt
    end
    nothing
end

function _get_step_vectors!(step_vectors::Vector{CartesianVector{T}}, Δt::T, cdm::AbstractChargeDriftModel{T}, ::Type{Hole})::Nothing where {T <: SSDFloat}
    for n in eachindex(step_vectors)
        step_vectors[n] = getVh(SVector{3,T}(step_vectors[n]), cdm) * Δt
    end
    nothing
end


function _check_and_update_position!(step_vectors::Vector{CartesianVector{T}}, 
            current_pos::Vector{CartesianPoint{T}},
            done::Vector{Bool},
            normal::Vector{Bool},
            drift_path::Array{CartesianPoint{T},2},
            timestamps::Vector{T},
            istep::Int,
            det::SolidStateDetector{T, S},
            g::Grid{T, 3, S},
            point_types::PointTypes{T, 3, S},
            startpos::Vector{CartesianPoint{T}},
            Δt::T,
            verbose::Bool
        )::Nothing where {T <: SSDFloat, S}

    for n in eachindex(normal)
        normal[n] = done[n] || _is_next_point_in_det(current_pos[n]+step_vectors[n], CylindricalPoint(current_pos[n]+step_vectors[n]), det, point_types)
    end
    
    if all(normal)
        #all charges are either finished or still inside the detector => drift normally
        current_pos .+= step_vectors
        drift_path[:,istep] .= current_pos
        timestamps[istep] = timestamps[istep-1] + Δt
    else
        #all charges that would not be inside after the drift step
        for n in findall(.!normal)
            crossing_pos::CartesianPoint{T}, cd_point_type::UInt8, boundary_index::Int, surface_normal::CartesianVector{T} = 
                get_crossing_pos(det, g, copy(current_pos[n]), current_pos[n] + step_vectors[n])
            if cd_point_type == CD_ELECTRODE
                done[n] = true
                drift_path[n,istep] = crossing_pos
                timestamps[istep] = timestamps[istep-1] + Δt      
            elseif cd_point_type == CD_FLOATING_BOUNDARY
                projected_vector::CartesianVector{T} = CartesianVector{T}(project_to_plane(step_vectors[n], surface_normal))
                projected_vector = modulate_surface_drift(projected_vector)
                next_pos::CartesianPoint{T} = current_pos[n] + projected_vector
                small_projected_vector = projected_vector * T(0.001)
                i::Int = 0
                while i < 1000 && !(next_pos in det)
                    next_pos -= small_projected_vector
                    i += 1
                end
                if i == 1000 && verbose @warn("Handling of charge at floating boundary did not work as intended. Start Position (Cart): $(startpos[n])") end
                drift_path[n,istep] = next_pos
                step_vectors *= (1 - i * T(0.001))  # scale down the step_vectors for all other charge clouds
                Δt *= (1 - i * T(0.001))            # scale down Δt for all charge clouds
                done[n] = is_zero_vector(geom_round.(next_pos - current_pos[n]))
                current_pos[n] = next_pos
            else # if cd_point_type == CD_BULK or CD_OUTSIDE
                if verbose @warn ("Internal error for charge starting at $(startpos[n])") end
                done[n] = true
                drift_path[n,istep] = current_pos[n]
                timestamps[istep] = timestamps[istep-1] + Δt
            end  
        end
        #drift all other charge clouds normally according to the new Δt_min
        for n in findall(normal)
            current_pos[n] += step_vectors[n]
            drift_path[n,istep] = current_pos[n]
        end 
    end
    nothing
end



"""
    _drift_charge!(...)

Before calling this function one should check that `startpos` is inside `det`: `in(startpos, det)`
"""
function _drift_charge!(
                            drift_path::Array{CartesianPoint{T},2},
                            timestamps::Vector{T},
                            det::SolidStateDetector{T, S},
                            point_types::PointTypes{T, 3, S},
                            g::Grid{T, 3, S},
                            startpos::Vector{CartesianPoint{T}},
                            charges::Vector{T},
                            Δt::T,
                            velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3},
                            cdm::AbstractChargeDriftModel{T},
                            CC::Type,
                            diffusion::Bool = false,
                            self_repulsion::Bool = false;
                            verbose::Bool = true
                        )::Int where {T <: SSDFloat, S}
    
    n_hits::Int = length(startpos)
    drift_path[:,1] = startpos
    timestamps[1] = T(0)
    max_nsteps::Int = size(drift_path,2)
    ϵr::T = T(det.semiconductors[1].material.ϵ_r)
    
    last_real_step_index::Int = 1
    current_pos::Vector{CartesianPoint{T}} = copy(startpos) # copy is needed to prevent overwriting startpos
    step_vectors::Vector{CartesianVector{T}} = Vector{CartesianVector{T}}(undef, n_hits)
    done::Vector{Bool} = fill(false, n_hits)
    normal::Vector{Bool} = fill(false, n_hits)
    
    for istep in 2:max_nsteps
        last_real_step_index += 1
        _set_to_zero_vector!(step_vectors)
        _add_fieldvector_drift!(step_vectors, current_pos, done, velocity_field, det)
        if diffusion _add_fieldvector_diffusion!(step_vectors) end
        if self_repulsion _add_fieldvector_selfrepulsion!(step_vectors, current_pos, done, charges, ϵr, CC) end
        _modulate_fieldvector!(step_vectors, current_pos, det.virtual_drift_volumes)    
        _get_step_vectors!(step_vectors, Δt, cdm, CC)
        _check_and_update_position!(step_vectors, current_pos, done, normal, drift_path, timestamps, istep, det, g, point_types, startpos, Δt, verbose)
        if all(done) break end
    end
    last_real_step_index
end

# Point types for charge drift: Defined in DetectorGeometries/DetectorGeometries.jl
# const CD_ELECTRODE = 0x00
# const CD_OUTSIDE = 0x01
# const CD_BULK = 0x02
# const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

function get_crossing_pos(  detector::SolidStateDetector{T, S}, grid::Grid{T, 3}, point_in::CartesianPoint{T}, point_out::CartesianPoint{T};
                            max_n_iter::Int = 500)::Tuple{CartesianPoint{T}, UInt8, Int, CartesianVector{T}} where {T <: SSDFloat, S}
    point_mid::CartesianPoint{T} = T(0.5) * (point_in + point_out)
    cd_point_type::UInt8, contact_idx::Int, surface_normal::CartesianVector{T} = point_type(detector, grid, _convert_vector(point_mid, Val(S)))
    for i in 1:max_n_iter
        if cd_point_type == CD_BULK
            point_in = point_mid
        elseif cd_point_type == CD_OUTSIDE
            point_out = point_mid
        elseif cd_point_type == CD_ELECTRODE
            break
        else #elseif cd_point_type == CD_FLOATING_BOUNDARY
            break
        end
        point_mid = T(0.5) * (point_in + point_out)
        cd_point_type, contact_idx, surface_normal = point_type(detector, grid, _convert_vector(point_mid, Val(S)))
    end
    return point_mid, cd_point_type, contact_idx, surface_normal
end


include("plot_recipes.jl")
