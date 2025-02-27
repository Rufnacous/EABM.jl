
# Type definitions such that we can have Articulation{mutable struct FluidborneFilamentProperties},
# allowing torque_elastic and fluid forces to be applied to an articulation.

abstract type FFIProperties end;

mutable struct StandardFFIProperties <: FFIProperties
    drag_coefficient::Number
    friction_coefficient::Number
    frontal_area::Number
    skin_area::Number
    volume::Number
    added_fluid_volume::Number
    added_inertia::Matrix{<: Number}
    fluid_freedom_subspace::Matrix{<: Number}
end


mutable struct FluidborneFilamentProperties
    elastics::ElasticProperties
    ffi::StandardFFIProperties
end