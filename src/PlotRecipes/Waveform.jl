"""
    get_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T}, contact_id::Int)
    
Returns the electron and hole contribution to the waveform of a [`Contact`](@ref) with a given
`contact_id` of an [`Event`](@ref) as a `NamedTuple` with two entries: 
`electron_contribution` and `hole_contribution`.

## Arguments
* `evt::Event{T}`: [`Event`](@ref) in which the charges have already been drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` were drifted.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the waveform should be split into electron and hole contribution.

## Example
```julia 
using Plots
using SolidStateDetector
T = Float32

simulation = Simulation{T}(SSD_examples[:InvertedCoax])
simulate!(simulation)
event = Event([CartesianPoint{T}(0.02,0.01,0.05)])
simulate!(event, simulation)

contact_id = 1
wf = get_electron_and_hole_contribution(evt, sim, contact_id)
```

!!! note 
    The drift paths in `evt` need to be calculated using [`drift_charges!`](@ref) before calling this function.
    
See also [`plot_electron_and_hole_contribution`](@ref).
"""
function get_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T, S}, contact_id::Int; waveform_unit::Unitful.Units = u"e_au",
            )::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} where {T <: SSDFloat, S}
    
    @assert !ismissing(evt.drift_paths) "The charge drift is not yet simulated. Please use `drift_charges!(evt, sim)`!"
    
    dt::T = T(ustrip(u"ns", diff(evt.drift_paths[1].timestamps_e)[1]*u"s"))
    wp::Interpolations.Extrapolation{T, 3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    signal_e::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_e, evt.drift_paths))))
    signal_h::Vector{T} = zeros(T, length(maximum(map(p -> p.timestamps_h, evt.drift_paths))))

    ctm = sim.detector.semiconductor.charge_trapping_model
    for i in eachindex(evt.drift_paths)
        energy = flatview(evt.energies)[i]
        
        dp_e::Vector{CartesianPoint{T}} = evt.drift_paths[i].e_path
        dp_e_t::Vector{T} = evt.drift_paths[i].timestamps_e
        add_signal!(signal_e, dp_e_t, dp_e, dp_e_t, -energy, wp, sim.point_types, ctm)
        
        dp_h::Vector{CartesianPoint{T}} = evt.drift_paths[i].h_path
        dp_h_t::Vector{T} = evt.drift_paths[i].timestamps_h
        add_signal!(signal_h, dp_h_t, dp_h, dp_h_t, energy, wp, sim.point_types, ctm)
    end
    calibration_factor::Quantity{T, dimension(waveform_unit)} = _convert_internal_energy_to_external_unit(waveform_unit, sim.detector.semiconductor.material)
    return (electron_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_e)), signal_e * calibration_factor),
            hole_contribution = RDWaveform(range(zero(T) * u"ns", step = dt * u"ns", length = length(signal_h)), signal_h * calibration_factor))
end

export get_electron_and_hole_contribution


"""
    plot_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T}, contact_id::Int)
    
Plots the waveform as well as the electron and hole contribution to the waveform 
of a [`Contact`](@ref) with a given `contact_id` of an [`Event`](@ref).

## Arguments
* `evt::Event{T}`: [`Event`](@ref) in which the charges have already been drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` were drifted.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the waveform should be split into electron and hole contribution.

## Keywords 
* `n_samples::Int`: Number of samples with which the waveforms will be plotted. The default is the number of samples of the original waveform.

## Example
```julia 
using Plots
using SolidStateDetector
T = Float32

simulation = Simulation{T}(SSD_examples[:InvertedCoax])
simulate!(simulation)
event = Event([CartesianPoint{T}(0.02,0.01,0.05)])
simulate!(event, simulation)

contact_id = 1
plot_electron_and_hole_contribution(evt, sim, contact_id, n_samples = 300)
```

!!! note 
    The drift paths in `evt` need to be calculated using [`drift_charges!`](@ref) before calling this function.
    
!!! note 
    This method requires to load the package [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
    
See also [`get_electron_and_hole_contribution`](@ref).
"""
function plot_electron_and_hole_contribution end

@userplot Plot_electron_and_hole_contribution
@recipe function f(gdd::Plot_electron_and_hole_contribution; linewidth = 2, n_samples = missing)
    
    sim::Simulation = gdd.args[2]
    T = get_precision_type(sim)
    evt::Event{T} = gdd.args[1]
    contact_id::Int = gdd.args[3]
    
    ismissing(n_samples) ? n_samples = length(evt.waveforms[contact_id].signal) : nothing

    # check in which units the waveforms are calibrated
    waveform_unit::Unitful.Units = unit(evt.waveforms[contact_id].signal[end])
    wf::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} = get_electron_and_hole_contribution(evt, sim, contact_id; waveform_unit)
    
    unitformat --> :slash
    yguide --> waveform_unit isa Unitful.Units{<:Any, Unitful.Charge} ? "Charge" : "Energy"

    @series begin
        linecolor := :red 
        linewidth --> linewidth
        label --> "Electron contribution"
        add_baseline_and_extend_tail(wf.electron_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor := :green 
        linewidth --> linewidth
        label --> "Hole contribution"
        add_baseline_and_extend_tail(wf.hole_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor --> contact_id
        linewidth --> linewidth
        seriesalpha --> 0.5
        label --> "Sum"
        add_baseline_and_extend_tail(evt.waveforms[contact_id], 0, n_samples)
    end
end


# This should be in RadiationDetectorSignals.jl
@recipe function f(wvs::Vector{Union{Missing, RadiationDetectorSignals.RDWaveform}})
    @series begin
        RadiationDetectorSignals.RDWaveform[wv for wv in skipmissing(wvs)]
    end
end