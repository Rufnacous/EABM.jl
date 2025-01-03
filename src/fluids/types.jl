
abstract type FFIProperties end;

mutable struct StandardFFIProperties <: FFIProperties
    drag_coefficient::Number
    frontal_area::Number
    added_fluid_volume::Number
    added_inertia::Matrix{<: Number}
    fluid_freedom_subspace::Matrix{<: Number}
end


mutable struct FluidborneFilamentProperties
    elastics::ElasticProperties
    ffi::StandardFFIProperties
end