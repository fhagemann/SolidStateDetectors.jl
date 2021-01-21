abstract type ParticleType end
abstract type Alpha <: ParticleType end
abstract type Beta <: ParticleType end
abstract type Gamma <: ParticleType end

radius_guess(charge::T, ::Type{Alpha}) where {T} = T(0.0001)
radius_guess(charge::T, ::Type{Beta}) where {T} = T(0.0005)
radius_guess(charge::T, ::Type{Gamma}) where {T} = T(0.0005) 
