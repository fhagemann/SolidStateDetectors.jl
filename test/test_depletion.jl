using SolidStateDetectors
using Test
using Unitful

T = Float32

@testset "Test depletion estimation" begin

    sim = Simulation{T}("BEGe_01.yaml")
    calculate_electric_potential!(sim, refinement_limits=0.01)
    id = SolidStateDetectors.determine_bias_voltage_contact_id(sim.detector)
    @test id == 1

    calculate_weighting_potential!(sim, id, refinement_limits=0.01)
    # SolidStateDetectors._adapt_weighting_potential_to_electric_potential_grid!(sim, id)
    U_est = estimate_depletion_voltage(sim) # around -2600V
    ΔU = 50u"V"

    # test slight underdepletion
    U₊ = U_est + ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₊)
    calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    @test !is_depleted(sim.point_types)
    
    # test slight overdepletion
    U₋ = U_est - ΔU
    sim.detector = SolidStateDetector(sim.detector, contact_id=id, contact_potential=U₋)
    calculate_electric_potential!(sim, refinement_limits=0.01, depletion_handling=true)
    @test is_depleted(sim.point_types)
end
